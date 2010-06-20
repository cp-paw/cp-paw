!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_YLMPOL(LX,YLMPOL)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX
      REAL(8)   ,INTENT(OUT):: YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2)
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     ==========================================================================
!     == R^L*Y_LM = SUM_IJK X^IY^JZ^K * YLMPOL(IJK,LM)                        ==
!     ==========================================================================
      YLMPOL(:,:)=0.D0
      YLMPOL(1,1)=1.D0/SQRT(4.D0*PI)
      IF(LX.EQ.0) RETURN
      YLMPOL(4,2)=SQRT(3.D0/(4.D0*PI))
      YLMPOL(2,3)=SQRT(3.D0/(4.D0*PI))
      YLMPOL(3,4)=SQRT(3.D0/(4.D0*PI))
      IF(LX.EQ.1) RETURN
      YLMPOL(10,5)=SQRT(15.D0/(16.D0*PI))
      YLMPOL(7,5)=-SQRT(15.D0/(16.D0*PI))
      YLMPOL(8,6)=SQRT(60.D0/(16.D0*PI))
      YLMPOL(5,7)=2.D0*SQRT(5.D0/(16.D0*PI))
      YLMPOL(7,7)=-SQRT(5.D0/(16.D0*PI))
      YLMPOL(10,7)=-SQRT(5.D0/(16.D0*PI))
      YLMPOL(6,8)=SQRT(60.D0/(16.D0*PI))
      YLMPOL(9,7)=SQRT(60.D0/(16.D0*PI))
      IF(LX.EQ.2) RETURN
      CALL ERROR$MSG('IMPLEMENTED ONLY UP TO LX=2')
      CALL ERROR$I4VAL('LX',LX)
      CALL ERROR$STOP('GAUSSIAN_YLMPOL')
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
      REAL(8)                  :: AX2(3),FAC(3),fac1(3),XJ
      INTEGER(4)               :: J,K,N,IND1,IND2,I1,I2,J1,J2,K1,K2,M1,M2,nmk
      REAL(8)                  :: GARR(0:NX,3)
      real(8)                  :: B(0:NX,0:NX) !BINOMIAL COEFFICIENTS
!     **************************************************************************
!     == without a shift, return the identity ==================================
      if(sum(r(:)**2).lt.1.d-10) then
        c(:,:)=0.d0
        do n=1,nijk
          c(n,n)=1.d0
        enddo
        return
      end if
!
!     ==========================================================================
!     == CALCULATE BINOMIAL COEFFICIENTS                                      ==
!     ==========================================================================
      B(:,:)=0.D0
      DO N=0,NX
        B(N,0)=1.D0
        DO K=1,N-1
          B(N,K)=B(N-1,K-1)+B(N-1,K)
        ENDDO
        B(N,N)=1.D0
      ENDDO
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
!     == EXPAND POLYNOMIAL ABOUT OTHER CENTER AND assemble RESULT             ==
!     ==========================================================================
      C1(:,:)=0.D0
      C2(:,:)=0.D0
      C3(:,:)=0.D0
      FAC(:)=1.D0
      DO K=0,NX
        DO N=K,NX
          FAC1(:)=FAC(:)*B(N,K)
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
!!$        do j=0,nx
!!$          k1=max(n+j-nx,0)
!!$          fac(:)=(-r(:))**k1
!!$          do k=k1,n
!!$            C1(n+j-k,n)=C1(n+J-K,n)+GARR(J,1)*B(N,K)*FAC(1)
!!$            C2(N+J-K,n)=C2(n+j-k,n)+GARR(J,2)*B(N,K)*FAC(2)
!!$            C3(N+J-K,n)=C3(n+j-k,n)+GARR(J,3)*B(N,K)*FAC(3) 
!!$            fac(:)=-r(:)*fac(:)
!!$          enddo
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
      INTEGER(4)              :: N,M,I,J,K,IND
      REAL(8)                 :: GX,GY,GZ
!     **************************************************************************
      IF(ID.EQ.'CARTESIAN') THEN
        g=0.d0
        do ind=1,nijk
          CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',ind,i,j,k)
          g=g+c(ind)*r(1)**i*r(2)**j*r(3)**k
        enddo
        g=g*exp(-e*sum(r(:)**2))
        return
      end if 
!
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJK,N,J,K)
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
      SUBROUTINE GAUSSIAN_GAUSSINDEX(ID,IND,I,J,K)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(INOUT) :: IND
      INTEGER(4)  ,INTENT(INOUT) :: I,J,K
      INTEGER(4)                 :: N,IND1,I1
!     **************************************************************************
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
      SUBROUTINE GAUSSIAN_FITGAUSS(GID,NR,W,L,F,NEP,NPOW,EP,C)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: W(NR) ! WEIGHTING FOR FIT
      REAL(8)   ,INTENT(IN) :: F(NR) ! FUNCTION TO BE FITTED
      INTEGER(4),INTENT(IN) :: NEP   ! #(EXPONENTS)
      INTEGER(4),INTENT(IN) :: NPOW  ! #(POWERS) HIGHEST POWER=L+2*(NPOW-1)
      REAL(8)   ,INTENT(IN) :: EP(NEP)
      REAL(8)   ,INTENT(OUT):: C(NPOW,NEP)
      REAL(8)               :: G(NR,NPOW,NEP)
      REAL(8)               :: R2(NR),RL(NR)
      REAL(8)               :: B(NPOW,NEP)
      REAL(8)               :: A(NPOW,NEP,NPOW,NEP)
      REAL(8)               :: svar,AUX(NR)
      REAL(8)               :: Q(2*NPOW,NEP,NEP)
      REAL(8)               :: scale(NPOW,NEP)
      REAL(8)               :: WR2(NR)
      INTEGER(4)            :: I,J,I1,J1,I2,J2
      logical(4),parameter  :: ttest=.true.
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R2)
      RL(:)=R2(:)**L
      R2(:)=R2(:)**2
      WR2(:)=W(:)*R2(:)
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
!     ==  CONSTRUCT B=<gi|W|Gj>                                               ==
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
      a(:,:,:,:)=0.d0
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
!     ==  scale test functions to avoid numerical problems
!     ==========================================================================
      do i=1,nep
        svar=exp(1.d0)/(2.d0*ep(i))
        do j=1,npow
          scale(j,i)=1.d0/(svar*real(l+2*j-2))**(l+2*j-2)
        enddo
      enddo
      DO I2=1,NEP
        DO J2=1,NPOW
          DO I1=1,NEP
            DO J1=1,NPOW
              A(J1,I1,J2,I2)=scale(j1,i1)*a(j1,i1,j2,i2)*scale(j2,i2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO I1=1,NEP
        DO J1=1,NPOW
          b(J1,I1)=scale(j1,i1)*b(j1,i1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  solve least squares equation                                        ==
!     ==========================================================================
      CALL LIB$MATRIXSOLVER8(NEP*NPOW,NEP*NPOW,1,A,C,B)
!
!     ==========================================================================
!     ==  scale coefficients
!     ==========================================================================
      DO I1=1,NEP
        DO J1=1,NPOW
          c(J1,I1)=scale(j1,i1)*c(j1,i1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  test
!     ==========================================================================
      if(ttest) then
        do i1=1,npow
          do j1=1,nep
            svar=0.d0
            do i2=1,npow
              do j2=1,nep
                svar=svar+a(i1,j1,i2,j2)*c(i2,j2)/scale(i2,j2)
              enddo
            enddo
            svar=svar-b(i1,j1)
print*,'svar ',l,i1,j1,svar
            if(abs(svar).gt.1.d-6) then
              print*,'fitting error ',svar,i1,j1
            end if
          enddo
        enddo

print*,'b/scale ',b/scale
print*,'scale ',scale
print*,'c/scale ',c/scale
print*,'c ',c
!
        AUX(:)=WR2(:)*F(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,svar)
        do i1=1,npow
          do j1=1,nep
            do i2=1,npow
              do j2=1,nep
                svar=svar+c(i1,j1)*a(i1,j1,i2,j2)*c(i2,j2) &
      &                  /scale(i1,j1)/scale(i2,j2)
              enddo
            enddo
            svar=svar-2.d0*b(i1,j1)/scale(i1,j1)*c(i1,j1)
          enddo
        enddo
print*,'square deviation ',svar
      end if
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
!     **************************************************************************
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKA',NIJKA)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKB',NIJKB)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKP,NP,J,K)
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
                  INDP=0
                  DO MP=0,NP
                    DO IP=0,MP
                      DO JP=0,MP-IP
                        KP=MP-IP-JP
                        INDP=INDP+1
                        CP(INDP)=CP(INDP)+CA(INDA)*CB(INDB) &
        &                       *HX(IA,IB,IP)*HY(JA,JB,JP)*HZ(KA,KB,KP)
                      ENDDO
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
      SUBROUTINE GAUSSIAN_OVERLAP(NIJKA,nca,EA,RA,CA,NIJKB,ncb,EB,RB,CB,SAB)
!     **************************************************************************
!     ** EVALUATES THE OVERLAP OF TWO SETS OF FUNCTIONS REPRESENTED BY        **
!     ** CARTESIAN GAUSSIANS.                                                 **
!     ** SET "A" WITH NCA FUNCTIONS HAS THE EXPONENT EA AND IS CENTERED AT RA **
!     ** SET "B" WITH NCB FUNCTIONS HAS THE EXPONENT EB AND IS CENTERED AT RB **
!     ** THE RESULTING OVERLAP MATRIX IS S_IJ=<A_I|B_J>                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: Nca
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: CA(NIJKA,nca)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: Ncb
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: CB(NIJKB,ncb)
      REAL(8)   ,INTENT(OUT) :: SAB(nca,ncb)
      REAL(8)               :: EP
      INTEGER(4)            :: NA,NB
      INTEGER(4)            :: NABX
      INTEGER(4)            :: i,J,k
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA
      INTEGER(4)            :: INDB,MB,IB,JB,KB
      real(8)               :: pi
      real(8)               :: svar
!     **************************************************************************
      pi=4.d0*atan(1.d0)
!
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
!
!     ==========================================================================
!     ==  determine hermite coefficients                                      ==
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
!     ==  work out overlap matrix                                             ==
!     ==========================================================================
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
                  svar=HX(IA,IB,0)*HY(JA,JB,0)*HZ(KA,KB,0)
                  do j=1,ncb
                    do i=1,nca
                      SAB(i,j)=SAB(i,j)+svar*CA(INDA,i)*CB(INDB,j)
                    enddo
                  enddo
!!$if(abs(svar*ca(inda,1)*cb(indb,1)/sab(1,1)).gt.1.d-1) then
!!$write(*,fmt='("sab(1,1) ",2i6,10e10.2)')inda,indb,sab(1,1),ca(inda,1),cb(indb,1)
!!$end if
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      SAB(:,:)=SAB(:,:)*(pi/(ea+eb))**1.5d0
!
!!$do i=1,nijka
!!$  write(*,fmt='("ca",i8,20e10.2)')i,ca(i,:),cb(i,:)
!!$enddo
!!$do j=1,ncb
!!$  write(*,fmt='("s",20e10.2)')sab(:,j)
!!$enddo
!!$stop
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEC(N,XA,XB,EA,EB,H)
!     **************************************************************************
!     ** EVALUATES A SET OF Hermite coefficients                              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: XA
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: XB
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(OUT):: H(0:N,0:N,0:2*N)
      REAL(8)               :: EP,MU,ONEBY2P,XP,PA,PB,TP1
      REAL(8)               :: PI
      INTEGER(4)            :: I,J,IJSUM,T,I1,I2
!     **************************************************************************
      IF(N.LT.1) THEN
        call error$STOP('GAUSSIAN_HERMITEC')
      END IF
      H(:,:,:)=0.D0
      EP=EA+EB
      MU=EA*EB/EP
      ONEBY2P=1.D0/(2.D0*EP)
      XP=(EA*XA+EB*XB)/EP
      PA=XP-XA
      PB=XP-XB
      H(0,0,0)=EXP(-MU*(XB-XA)**2)
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
      SUBROUTINE GAUSSIAN_FOURCENTER(NIJKA,EA,RA,CA,NIJKB,EB,RB,CB &
     &                              ,NIJKC,EC,RC,CC,NIJKD,ED,RD,CD,UABCD)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NIJKC
      INTEGER(4),INTENT(IN) :: NIJKD
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(IN) :: EC
      REAL(8)   ,INTENT(IN) :: ED
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: RC(3)
      REAL(8)   ,INTENT(IN) :: RD(3)
      REAL(8)   ,INTENT(IN) :: CA(3)
      REAL(8)   ,INTENT(IN) :: CB(3)
      REAL(8)   ,INTENT(IN) :: CC(3)
      REAL(8)   ,INTENT(IN) :: CD(3)
      REAL(8)   ,INTENT(OUT) :: UABCD
      INTEGER(4)             :: NIJKP,NIJKQ
      REAL(8)                :: EP,EQ
      REAL(8)                :: RP(3),RQ(3)
      REAL(8)   ,ALLOCATABLE :: CP(:),CQ(:)
      INTEGER(4)             :: NA,NB,NC,ND,NP,NQ,N
      INTEGER(4)             :: J,K
      REAL(8)   ,ALLOCATABLE :: HI(:,:,:)
      REAL(8)                :: SIGN
      INTEGER(4)             :: INDP,INDQ,MP,MQ,IP,IQ,JP,JQ,KP,KQ
      REAL(8)                :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSIAN_FOURCENTER')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSIAN_FOURCENTER')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKC,NC,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSIAN_FOURCENTER')
      END IF
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKD,ND,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSIAN_FOURCENTER')
      END IF
!
      NP=NA+NB 
      NQ=NC+ND 
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKQ,NQ,0,0)
      ALLOCATE(CP(NIJKP))
      ALLOCATE(CQ(NIJKQ))
      CALL GAUSSIAN_CONTRACT2GAUSS(NIJKA,EA,RA,CA,NIJKB,EB,RB,CB &
     &                            ,NIJKP,EP,RP,CP)
      CALL GAUSSIAN_CONTRACT2GAUSS(NIJKC,EC,RC,CC,NIJKD,ED,RD,CD &
     &                            ,NIJKQ,EQ,RQ,CQ)
      N=MAX(NP,NQ)
      ALLOCATE(HI(0:N,0:N,0:N))
      CALL GAUSSIAN_HERMITEINTEGRAL(N,EP*EQ/(EP+EQ),RP-RQ,HI)
      UABCD=0.D0
      INDP=0
      DO MP=0,NP      
        DO IP=0,MP
          DO JP=0,MP-IP
            KP=MP-IP-JP
            INDP=INDP+1
!         
            INDQ=0
            DO MQ=0,NQ      
              DO IQ=0,MQ
                DO JQ=0,MQ-IQ
                  KQ=MQ-IQ-JQ
                  INDQ=INDQ+1
                  SIGN=IQ+JQ+KQ
                  UABCD=UABCD+SIGN*HI(IP+IQ,JP+JQ,KP+KQ) &
     &                 *CP(INDP)*CQ(INDQ)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      UABCD=(2.D0*PI)**2.5D0/(EP*EQ*SQRT(EP+EQ))*UABCD
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEINTEGRAL(N,E,R,HI)
!     **************************************************************************
!     ** EVALUATES A HERMITE COULOMB INTEGRALS                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(OUT):: HI(0:N,0:N,0:N)
      REAL(8)               :: F(0:N)
      INTEGER(4)            :: T,U,V
      INTEGER(4)            :: M
      REAL(8)               :: X
!     **************************************************************************
!
!     ==========================================================================
!     ==  DETERMINE BOYS' FUNCTIONS                                           ==
!     ==========================================================================
      X=E*(R(1)**2+R(2)**2+R(3)**2)
      CALL GAUSSIAN_BOYS(N,X,F)
      DO M=1,N
        F(M)=(-2.D0*E)*F(M)    !R^N_{000}
      ENDDO
!
!     ==========================================================================
!     ==  USE RECURSION TO CONSTRUCT HERMITE COULOMB INTEGRALS                ==
!     ==========================================================================
      DO M=1,3*N
        HI(0,0,0)=F(M)
        DO T=MIN(M,N),2,-1
          DO U=0,M-T
            V=M-T-U
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
!     **                                                                      **
!     ** METHOD USED:                                                         **
!     **   B.A. MAMEDOV, "ON THE EVALUATION OF BOYS FUNCTIONS USING DOWNWARD  **
!     **   RECURSION RELATION", J. MATH. CHEM. 36, P301                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: X
      REAL(8)   ,INTENT(OUT):: F(0:N)
      REAL(8)   ,PARAMETER  :: D=20.D0   ! #(ACCURATE DIGITS) 
      REAL(8)               :: Y,EXPMX,PI
      INTEGER(4)            :: MT,M
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      IF(ABS(X-REAL(N)).GT.1.D-5) THEN
        MT=N+D/ABS(LOG10(REAL(N)/X))
      ELSE
        MT=N+D/ABS(LOG10(REAL(N)))
      END IF   
      MT=2*INT(0.5D0*REAL(MT+2))
      EXPMX=EXP(-X)
      Y=SQRT(0.25*PI/X)
      DO M=MT,0,-1
        Y=(2.D0*X*Y+EXPMX)/REAL(2*M+1,KIND=8) !=F(M-1)
        IF(M.LE.N)F(M)=Y
      ENDDO
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
      INTEGER(4),PARAMETER :: NFIL=10
      INTEGER(4),PARAMETER :: N=50
      REAL(8)              :: F(0:N)
      REAL(8)              :: X
      INTEGER(4)           :: M
      REAL(8)              :: XA,XB,EA,EB
      REAL(8)              :: H(0:N,0:N,0:2*N)
      INTEGER(4)           :: I,J
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
      SUBROUTINE GAUSSIAN_TEST_FOURCENTER()
!     **************************************************************************
!     ** TESTROUTINE FOR GAUSSIAN_FOURCENTER                                  **
!     ** CURRENTLY IT JUST CALCULATES A MATRIX ELEMENT BETWEEN FUNCTIONS      **
!     ** EXPRESSED AS EXPANSIONS IN CARTESIAN GAUSSIANS                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: NIJKA,NIJKB,NIJKC,NIJKD,NIJKP,NIJKQ
      REAL(8)    :: EA,EB,EC,ED,EP,EQ
      REAL(8)    :: RA(3),RB(3),RC(3),RD(3),RP(3),RQ(3)
      REAL(8),ALLOCATABLE:: CA(:),CB(:),CC(:),CD(:),CP(:),CQ(:)
      REAL(8)    :: UABCD
      INTEGER(4) :: NA,NB,NC,ND,NP,NQ
      INTEGER(4) :: I
!     **************************************************************************
      NA=2
      NB=2
      NC=2
      ND=2
      NP=9  ! RESULTS BECOME VERY POOR IF NP<NA+NB!!
      NQ=9  ! RESULTS BECOME VERY POOR IF NQ<NC+ND!!
      CALL RANDOM_NUMBER(EA)
      CALL RANDOM_NUMBER(EB)
      CALL RANDOM_NUMBER(EC)
      CALL RANDOM_NUMBER(ED)
      DO I=1,3
        CALL RANDOM_NUMBER(RA(I))
        CALL RANDOM_NUMBER(RB(I))
        CALL RANDOM_NUMBER(RC(I))
        CALL RANDOM_NUMBER(RD(I))
      ENDDO
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKA,NA,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKB,NB,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKC,NC,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKD,ND,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKQ,NQ,0,0)
      ALLOCATE(CA(NIJKA))
      ALLOCATE(CB(NIJKB))
      ALLOCATE(CC(NIJKC))
      ALLOCATE(CD(NIJKD))
      ALLOCATE(CP(NIJKP))
      ALLOCATE(CQ(NIJKQ))
      DO I=1,NIJKA
        CALL RANDOM_NUMBER(CA(I))
      ENDDO
      DO I=1,NIJKB
        CALL RANDOM_NUMBER(CB(I))
      ENDDO
      DO I=1,NIJKC
        CALL RANDOM_NUMBER(CC(I))
      ENDDO
      DO I=1,NIJKD
        CALL RANDOM_NUMBER(CD(I))
      ENDDO
 
      CALL GAUSSIAN_FOURCENTER(NIJKA,EA,RA,CA,NIJKB,EB,RB,CB &
     &                        ,NIJKC,EC,RC,CC,NIJKD,ED,RD,CD,UABCD)
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
      INTEGER(4)           :: I,IPROBE
      REAL(8)              :: RAN
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
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKA,NA,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKB,NB,0,0)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,0,0)
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
      INTEGER(4)           :: N
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
      INTEGER(4),PARAMETER :: NFIL=10
      INTEGER(4),PARAMETER :: N1=12,N2=6,N3=6
      REAL(8)              :: F1(N1),F2(N2),F3(N3)
      REAL(8)              :: X1(N1),X2(N2),X3(N3)
      INTEGER(4)           :: M1(N1),M2(N2),M3(N3)
      REAL(8)              :: F
      REAL(8)              :: DEV
      REAL(8)   ,PARAMETER :: TOL=1.D-8
      INTEGER(4)           :: I
!     **************************************************************************
!     == TABLE 1 ===============================================================
      M1=(/8,15,20,25,31,11,42,75,100,20,45,100/)
      X1=(/16.D0,27.D0,30.D0,13.D0,34.D0,38.D0,32.D0,30.D0,33.D0 &
     &   ,1.4D-3,6.4D-5,2.6D-7/)
      F1=(/4.02308592502660D-07,1.08359515555596D-11,1.37585444267909D-03 &
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
        DEV=ABS((F-F1(I))/F1(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M1(I),X1(I),F,F1(I),DEV
        IF(DEV.GT.TOL) THEN
          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
!
!     ==========================================================================
      DO I=1,N2
        CALL GAUSSIAN_BOYS(M2(I),X2(I),F)
        DEV=ABS((F-F2(I))/F2(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M2(I),X2(I),F,F2(I),DEV
        IF(DEV.GT.TOL) THEN
          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
!
!     ==========================================================================
      DO I=1,N3
        CALL GAUSSIAN_BOYS(M3(I),X3(I),F)
        DEV=ABS((F-F3(I))/F3(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M3(I),X3(I),F,F3(I),DEV
        IF(DEV.GT.TOL) THEN
          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
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
      REAL(8)              :: dR(3)
      REAL(8)              :: R(3)
      REAL(8)              :: T(NIJK,NIJK)
      REAL(8)              :: C0(NIJK)
      REAL(8)              :: C1(NIJK)
      REAL(8)              :: F1,F0
      INTEGER(4)           :: NFIL=6
      INTEGER(4)           :: I
!     **************************************************************************
      E=1.D0
      dR(:)=(/1.D0,0.D0,0.D0/)
      C0(:)=0.D0
      C0(4)=1.D0
!
      CALL GAUSSIAN_SHIFTCENTER(NX,NIJK,E,dR,T)
      C1=MATMUL(T,C0)
!
      open(nfil,file='dat',form='formatted')
      DO I=1,NP
        R(:)=0.D0
        R(1)=-3.D0+6.D0*REAL(I-1)/REAL(NP-1)
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E,C0,R-dr,F0)
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E,C1,R,F1)
        WRITE(NFIL,*)R(1),F0,F1,F1-F0
      ENDDO
      close(nfil)
      STOP
      RETURN
      END
