!***********************************************************************
!**                                                                   **
!**  NAME: SPHERICAL                                                  **
!**                                                                   **
!**  PURPOSE: REAL SPHERICAL HARMONICS AND CLEBSCH GORDON COEFFICIENTS**
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    GETYLM                                                         **
!**    ROTATEYLM                                                      **
!**    CLEBSCH                                                        **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) FUNCTIONS NEED TO BE PUT INTO THE FORM SPHERICAL$FUNCTION   **
!**    2) INTERNAL ROUTINES SHALL BE INCLUDED INTO A MODULE           **
!**                                                                   **
!***********************************************************************
!
!     .....................................................GETYLM.......
      SUBROUTINE GETYLM(LMX,R,YLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE REAL SPHERICAL HARMONICS YLM AT POINT R           **
!     **  THE SPHERICAL HARMONICS ARE ORDERD WITH INDICES             **
!     **        LM(L,M)=1+L**2+L-M                                    **
!     **                                                              **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:                **
!     **      YLM(1)=SQRT( 1/( 4*PI))    * 1                          **
!     **      YLM(2)=SQRT( 3/( 4*PI))    * X / R                      **
!     **      YLM(3)=SQRT( 3/( 4*PI))    * Z / R                      **
!     **      YLM(4)=SQRT( 3/( 4*PI))    * Y / R                      **
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2      **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2      **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2      **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2      **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2      **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LMX
      REAL(8),   INTENT(IN)   :: R(3)
      REAL(8),   INTENT(OUT)  :: YLM(LMX)
      COMPLEX(8)              :: EIPHI,EIMPHI
      REAL(8)                 :: PI,FPI,SQ2
      INTEGER(4)              :: LX
      REAL(8)                 :: DIS,DISXY
      REAL(8)                 :: COSTHE,SINPHI,COSPHI
      INTEGER(4)              :: L,LM0,M,LM
      REAL(8)                 :: FAC,SVAR
      INTEGER(4)              :: LMM,LMP
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FPI=4.D0*PI
      SQ2=DSQRT(2.D0)
      IF(LMX.LE.1) THEN
        IF(LMX.LE.0) RETURN
        YLM(1)=1.D0/DSQRT(FPI)
        RETURN
      END IF
      LX=INT(DSQRT(DBLE(LMX-1)+0.01D0))
      DISXY=DSQRT(R(1)**2+R(2)**2)
      DIS=DSQRT(DISXY**2+R(3)**2)
      IF(DABS(DIS).LT.1.D-12) THEN
        YLM(1)=1.D0/DSQRT(FPI)
        DO LM=2,LMX
          YLM(LM)=0.D0
        ENDDO
        RETURN
      ENDIF
      COSTHE=R(3)/DIS
      CALL PLGNDR(LMX,LX,COSTHE,YLM)
      IF(DISXY.NE.0.D0) THEN
        SINPHI=R(2)/DISXY
        COSPHI=R(1)/DISXY
        EIPHI=CMPLX(COSPHI,SINPHI,8)
      ELSE
        EIPHI=(1.D0,0.D0)
      END IF
      DO L=0,LX
        LM0=L*(L+1)+1
        FAC=DSQRT(DBLE(2*L+1)/FPI)
        YLM(LM0)=YLM(LM0)*FAC
        EIMPHI=(1.D0,0.D0)
        DO M=1,L
          EIMPHI=EIMPHI*EIPHI
          LMM=LM0+M
          LMP=LM0-M
          SVAR=FAC*YLM(LMM)*SQ2
          YLM(LMP)=SVAR*REAL(EIMPHI,KIND=8)
          YLM(LMM)=SVAR*AIMAG(EIMPHI)
        ENDDO
      ENDDO
      RETURN
      END
!
!     .....................................................PLGNDR.......
      SUBROUTINE PLGNDR(LMX,LX,X,PLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE THE ORDINARY LEGENDRE POLYNOMIALS                 **
!     **  TIMES THE FACTOR DSQRT((L-M)!/(L+M)!)                       **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: LMX
      INTEGER(4)  ,INTENT(IN)  :: LX
      REAL(8)     ,INTENT(IN)  :: X
      REAL(8)     ,INTENT(OUT) :: PLM(LMX)
      REAL(8)                  :: FACT,FAC2
      REAL(8)                  :: SVAR
      INTEGER(4)               :: LM,L,M,I1,LMM,LMP
      REAL(8)                  :: PMM,CM,C0,CP
!     ******************************************************************
      FACT=DSQRT(1.D0-X**2)
      IF(FACT.EQ.0.D0) THEN
        SVAR=-X
        LM=0
        DO L=0,LX
          SVAR=-SVAR
          DO M=1,2*L+1
            LM=LM+1
            PLM(LM)=SVAR
          ENDDO
        ENDDO
      END IF
      PMM=1.D0
      DO M=0,LX
        CM=0.D0
        C0=PMM
        LM=M**2+1
        PLM(LM)=C0
        DO L=M+1,LX
          LM=LM+2*L
          CP=(X*DBLE(2*L-1)*C0-DBLE(L+M-1)*CM)/DBLE(L-M)
          PLM(LM)=CP
          CM=C0
          C0=CP
        ENDDO
        PMM=-DBLE(2*M+1)*FACT*PMM
      ENDDO
!
      DO L=1,LX
        I1=L*(L+1)+1
        FACT=1.D0
        FAC2=1.D0
        DO M=1,L
          LMM=I1+M
          LMP=I1-M
          FACT=FACT*DSQRT(DBLE((L-M+1)*(L+M)))
          FAC2=-FAC2
          PLM(LMP)=PLM(LMP)/FACT
          PLM(LMM)=FAC2*PLM(LMP)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ROTATEYLM(LMX,ROT,YLMROT)
!     **                                                              **
!     **  PRODUCES A LMX*LMX MATRIX YLMROT                            **
!     **  THAT TRANSFORMS A COEFFICIENT VECTOR C_LM                   **
!     **  OF A FUNCTION EXPRESSED AS F(R)=SUM_LM Y_LM(R)*C_LM         **
!     **  INTO A COORDINATE SYSTEM R'=ROT*R                           **
!     **  WITH Y'_LM(R')=Y_LM(R)                                      **
!     **  SUCH THAT    F(R)=SUM_LM Y'_LM(R)*C'_LM                     **
!     **  WITH C'_LM1=SUM_LM2 YLMROT_LM1,LM2 C_LM2                    **
!     **  REMAINS INVARIANT UNDER ROTATION                            **
!     **                                                              **
!     **  WORKS ONLY FOR REAL HARMONICS  WITH Y_2=C*X;Y_3=C*Z;Y_4=C*Y **
!     **  USES THE CLEBSCH GORDAN ROUTINE CLEBSCH                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX             ! #(L,M-ANGULAR MOMENTA)
      REAL(8)   ,INTENT(IN) :: ROT(3,3)        ! ROTATION MATRIX FOR POSITIONS
      REAL(8)   ,INTENT(OUT):: YLMROT(LMX,LMX) ! ROTATION MATRIX FOR COEFFICIENTS
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.     ! PRINTS RESULT
      LOGICAL(4),PARAMETER  :: TTEST=.TRUE.    ! CHECKS WHETHER RO IS UNITARY
      INTEGER(4)            :: LX
      INTEGER(4)            :: I1,I2,J1,J2
      INTEGER(4)            :: INDEX(3)
      INTEGER(4)            :: LM1A,LM1B,LM2A,LM2B
      INTEGER(4)            :: LM1,LM2,LM3,LM4,LM5,LM6,L
      REAL(8)               :: SVAR,SVAR1,SVAR2,CG
!     ******************************************************************
!
!      =================================================================
!      == TEST WHETHER ROTATION MATRIX IS UNITARY                     ==
!      =================================================================
       IF(TTEST) THEN
         SVAR=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(2,3)*ROT(3,2)) &
      &      +ROT(1,2)*(ROT(2,3)*ROT(3,1)-ROT(2,1)*ROT(3,3)) &
      &      +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))
         IF(DABS(SVAR-1.D0).GT.1.D-12) THEN
           CALL ERROR$MSG('ROTATION MATRIX NOT UNITARY')
           CALL ERROR$STOP('YLMROT')
         END IF
       END IF
!
!      =================================================================
!      == INITIALIZE YLMROT FOR L=0 AND L=1                           ==
!      =================================================================
       LX=INT(SQRT(REAL(LMX)))-1
       YLMROT(:,:)=(0.D0,0.D0)
       IF(LX.GE.0)YLMROT(1,1)=1.D0
       IF(LX.GE.1) THEN
         INDEX(1)=2
         INDEX(2)=4
         INDEX(3)=3
         DO I1=1,3
           I2=INDEX(I1)
           DO J1=1,3
             J2=INDEX(J1)
             YLMROT(I2,J2)=ROT(I1,J1)
           ENDDO
         ENDDO
       END IF
!
!      =================================================================
!      == APPLY RECURSION FOR YLMROT FOR L>1                          ==
!      =================================================================
       DO L=2,LX
         LM1A=L**2+1
         LM1B=(L+1)**2
         LM2A=(L-1)**2+1
         LM2B=L**2
         SVAR2=0.D0
         DO LM3=2,4
           DO LM5=LM2A,LM2B
             CALL CLEBSCH(LM1A,LM3,LM5,CG)
             SVAR2=SVAR2+CG**2
           ENDDO
         ENDDO
         
         DO LM1=LM1A,LM1B
           DO LM2=LM1A,LM1B
             SVAR=0.D0
             DO LM3=2,4
               DO LM5=LM2A,LM2B
                 SVAR1=0.D0
                 DO LM4=2,4
                   DO LM6=LM2A,LM2B
                     CALL CLEBSCH(LM2,LM4,LM6,CG)
                     SVAR1=SVAR1+YLMROT(LM3,LM4)*YLMROT(LM5,LM6)*CG
                   ENDDO
                 ENDDO
                 CALL CLEBSCH(LM1,LM3,LM5,CG)
                 SVAR=SVAR+SVAR1*CG
               ENDDO
             ENDDO
             YLMROT(LM1,LM2)=SVAR/SVAR2
           ENDDO
         ENDDO 
       ENDDO
!
!      =================================================================
!      == PRINT RESULT FOR TEST PURPOSES                              ==
!      =================================================================
       IF(TPR) THEN
         DO L=0,LX
           LM1A=L**2+1
           LM1B=(L+1)**2
           DO LM1=LM1A,LM1B
             WRITE(*,FMT='(9F10.3)')YLMROT(LM1A:LM1B,LM1)
           ENDDO
         ENDDO
       END IF
       RETURN
       END
!
!     .......................................................................
      subroutine spherical$ylmtrans(lmx,c)
!     **                                                                   **
!     **  transformation matrix from real spherical harmonics Ybar         **
!     **  to angular momentum eigenstates y                                **
!     **                                                                   **
!     **   1=\SUM_{LM,LM'} |y_{LM}><y_{LM}|YBAR_{LM'}><yBAR_{LM'}|         **
!     **    =\SUM_{LM,LM'} |y_{LM}>C_{LM,LM'}<yBAR_{LM'}|                  **
!     **                                                                   **
      implicit none
      integer(4),intent(in) :: lmx
      complex(8),intent(out):: c(lmx,lmx)
      integer(4)            :: l,m,lm,lmp,lmm,I
      integer(4)            :: lmax
      real(8)               :: fac
      real(8)               :: sqrtinv
      LOGICAL(4)            :: TTEST=.FALSE.
      COMPLEX(8)            :: CTEST(LMX,LMX)
      REAL(8)               :: SVAR
!     ***********************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((Lmax+1)**2.NE.LMX) THEN
        call error$MSG('LMX MUST SPAN FULL L-SHELLS')
        call error$i4VAL('LMX',LMX)
        CALL ERROR$STOP('SPHERICAL_ylmtrans')
      END IF
      sqrtinv=1.d0/sqrt(2.d0)
      c=(0.d0,0.d0)
      do L=0,LMAX
        lmp=l**2+l+1
        lmm=lmp
        c(lmp,lmp)=(1.d0,0.D0)
        fac=1.d0
        do m=1,l
          fac=-fac
          lmp=lmp+1
          lmm=lmm-1
          c(lmp,lmp)=    (1.d0, 0.d0)*sqrtinv
          c(lmp,lmm)=    (0.d0,-1.d0)*sqrtinv
          c(lmm,lmp)=fac*(1.d0, 0.d0)*sqrtinv
          c(lmm,lmm)=fac*(0.d0, 1.d0)*sqrtinv
        enddo
      enddo
!
      IF(TTEST) THEN
        CTEST=MATMUL(TRANSPOSE(CONJG(C)),C)
        DO I=1,LMX
          CTEST(I,I)=CTEST(I,I)-(1.D0,0.D0)
        ENDDO
        SVAR=MAXVAL(ABS(CTEST))
        IF(SVAR.GT.1.D-10) THEN
          CALL ERROR$MSG('TRANSFORMATION IS NOT UNITARY')   
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('SPHERICAL_YLMTRANS')   
        END IF
      END IF
      return
      end
!
!     .......................................................................
      subroutine spherical$ER(lmx,x,y,z)
!     **                                                                   **
!     **  matrix elements of the unit vectors x/r, y/r, z/y                **
!     **                                                                   **
!     **       xi/r=\sum_{lm,lm'} |y_lm><y_lm|xi/r|y_lm'><ybar_lm'|        **
!     **           =\sum_{lm,lm'} |y_lm> xi_{lm,lm'} <ybar_lm'|            **
!     **                                                                   **
!     **       x_{lm,lm'}=<lm|x/r|lm'>                                     **
!     **       y_{lm,lm'}=<lm|y/r|lm'>                                     **
!     **       z_{lm,lm'}=<lm|z/r|lm'>                                     **
!     **                                                                   **
      implicit none
      integer(4),intent(in) :: lmx
      real(8)   ,intent(out):: x(lmx,lmx)
      real(8)   ,intent(out):: y(lmx,lmx)
      real(8)   ,intent(out):: z(lmx,lmx)
      real(8)               :: pi,sq4piby3
      integer(4)            :: lm1,lm2
!     ***********************************************************************
      pi=4.d0*datan(1.d0)
      sq4piby3=sqrt(4.d0*pi/3.d0)
      DO LM1=1,LMX
        DO LM2=1,LMX
          CALL CLeBSCH(LM1,LM2,2,x(lm1,lm2))
          CALL CLeBSCH(LM1,LM2,3,y(LM1,LM2))
          CALL CLeBSCH(LM1,LM2,4,z(LM1,LM2))
        ENDDO
      ENDDO
      x=x*SQ4PIBY3
      y=y*SQ4PIBY3
      z=z*SQ4PIBY3
      RETURN
      END
!
!     .......................................................................
      subroutine spherical$L(lmx,lx,ly,lz)
!     **                                                                   **
!     **  angular momentum matrix elements in a representation             **
!     **  of real spherical harmonics                                      **
!     **                                                                   **
!     **    lx_{lm,lm'}=<lm|lx|lm'>                                        **
!     **                                                                   **
      implicit none
      integer(4),intent(in) :: lmx
      complex(8),intent(out):: lX(lmx,lmx)
      complex(8),intent(out):: lY(lmx,lmx)
      complex(8),intent(out):: lz(lmx,lmx)
      complex(8)            :: c(lmx,lmx)
      INTEGER(4)            :: Lmax
      INTEGER(4)            :: LM,L,M
      INTEGER(4)            :: LOX(LMX),MOX(LMX)
      REAL(8)               :: LPFAC(LMX),LMFAC(LMX),LZFAC(LMX)
      COMPLEX(8)            :: LPMAT(LMX,LMX),LMMAT(LMX,LMX)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
!     ***********************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((Lmax+1)**2.NE.LMX) THEN
        CALL ERROR$STOP('SPHERICAL_L')
      END IF
      LM=0
      DO L=0,LMAX
        DO M=-L,L
          LM=LM+1
          LOX(LM)=L
          MOX(LM)=M
          LPFAC(LM)=SQRT(REAL((L-M)*(L+M+1),KIND=8))
          LMFAC(LM)=SQRT(REAL((L+M)*(L-M+1),KIND=8))
          LZFAC(LM)=REAL(M,KIND=8)
        ENDDO
      ENDDO
!
!     =======================================================================
!     == angular momenta in angular momentum eigenstates                   ==
!     =======================================================================
      lz=(0.d0,0.d0)
      lpMAT=(0.d0,0.d0)
      lmMAT=(0.d0,0.d0)
      do lm=1,lmx
        lz(lm,lm)=cmplx(lzfac(lm),0.d0)
        if(mox(lm).ne.+lox(lm))lpMAT(lm+1,lm)=CMPLX(lpfac(lm),0.D0)
        if(mox(lm).ne.-lox(lm))lmMAT(lm-1,lm)=CMPLX(lmfac(lm),0.D0)
      enddo
!
!     =======================================================================
!     == transform from a basis of angular momentum eigenstates to a basis ==
!     == of real spherical harmonics                                       ==
!     =======================================================================
      call spherical$ylmtrans(lmx,c)
      lz=matmul(lz,c)
      lpmat=matmul(lpmat,c)
      lmmat=matmul(lmmat,c)
      c=transpose(conjg(c))
      lz=matmul(c,lz)
      lpmat=matmul(c,lpmat)
      lmmat=matmul(c,lmmat)
!
!     =======================================================================
!     == cartesian angular momenta                                         ==
!     =======================================================================
      lx= 0.5d0   *(Lpmat+lmmat)
      lY=-0.5d0*CI*(Lpmat-lmmat)
      return
      end
!
!     .......................................................................
      subroutine SPHERICAL_TEST()
!     **                                                                   **
      implicit none
      INTEGER(4),PARAMETER :: LMX=9
      COMPLEX(8)           :: LX(LMX,LMX),LY(LMX,LMX),LZ(LMX,LMX)
      complex(8)           :: c(lmx,lmx)
      complex(8)           :: cdagger(lmx,lmx)
      complex(8)           :: mat(lmx,lmx)
      integer(4)           :: i,LM,L,M
      integer(4)           :: LMAX
      REAL(8)              :: SVAR
      integer(4)           :: LOX(LMX),MOX(LMX)
      complex(8)           :: ci=(0.d0,1.d0)
!     ***********************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((Lmax+1)**2.NE.LMX) THEN
        CALL ERROR$STOP('SPHERICAL_L')
      END IF
      LM=0
      DO L=0,LMAX
        DO M=-L,L
          LM=LM+1
          LOX(LM)=L
          MOX(LM)=M
         ENDDO
      ENDDO
!
      CALL spherical$L(lmx,lx,ly,lz)
      call spherical$ylmtrans(lmx,c)
      cdagger=transpose(conjg(c))
!
!     == CHECK EIGENVALUE EQUATION FOR LZ
      MAT=MATMUL(LZ,CDAGGER)
      DO I=1,LMX
        MAT(:,I)=MAT(:,I)-CDAGGER(:,I)*REAL(MOX(I))
      ENDDO
      SVAR=MAXVAL(ABS(MAT))
      PRINT*,'DEVIATION EIGENVALUE LZ ',SVAR

!     == CHECK EIGENVALUE EQUATION FOR L**2
      mat=matmul(lx,lx)+matmul(ly,ly)+matmul(lz,lz)
      PRINT*,'L**2'
      do i=1,lmx
        write(*,fmt='(20("(",2f6.2,")"))')mat(I,:)
      enddo
!
!     == PRINT TRANSFORMATION FROM REAL HARMONICS
      print*,'c'
      mat=c
      do i=1,lmx
        write(*,fmt='(20("(",2f6.2,")"))')mat(I,:)
      enddo
!
!     == CHECK IF C IS UNITARY
      print*,'cdagger*c'
      mat=matmul(cdagger,c)
      do i=1,lmx
        write(*,fmt='(20("(",2f6.2,")"))')mat(:,i)
      enddo
      print*,'c*cdagger'
      mat=matmul(c,cdagger)
      do i=1,lmx
        write(*,fmt='(20("(",2f6.2,")"))')mat(:,i)
      enddo
!
      print*,'check commutator relations'
      mat=matmul(lx,ly)-matmul(ly,lx)-ci*lz
      print*,'[lx,ly]-i*lz',maxval(abs(mat))
      mat=matmul(ly,lz)-matmul(lz,ly)-ci*lx
      print*,'[ly,lz]-i*lx',maxval(abs(mat))
      mat=matmul(lz,lx)-matmul(lx,lz)-ci*ly
      print*,'[lz,lx]-i*ly',maxval(abs(mat))
      RETURN
      END
!
!..........................................................CLBSCH........
MODULE CLEBSCH_MODULE
INTEGER(4)          :: LMXX=0
REAL(8),ALLOCATABLE :: CGMAT(:,:,:)
END MODULE CLEBSCH_MODULE
!
!     ....................................................CLBSCH........
      SUBROUTINE CLEBSCH(LM1,LM2,LM3,CG)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE CLEBSCH GORDAN COEFFICIENTS FOR REAL SPHERICAL    **
!     **  HARMONICS.                                                  **
!     **                                                              **
!     **    Y(LM1)*Y(LM2)= SUM(LM3|CG(LM1,LM2,LM3)*Y(LM3))            **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      USE CLEBSCH_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1
      INTEGER(4),INTENT(IN) :: LM2
      INTEGER(4),INTENT(IN) :: LM3
      REAL(8)   ,INTENT(OUT):: CG
      INTEGER(4)            :: LM123X
      INTEGER(4)            :: LMAX
!     ******************************************************************
      LM123X=MAX(LM1,LM2,LM3)
      IF(LM123X.GT.LMXX) THEN
        LMAX=INT(SQRT(REAL(LM123X-1)+0.01))
        LMAX=MAX(4,LMAX)
        LMXX=(LMAX+1)**2
        IF(ALLOCATED(CGMAT))DEALLOCATE(CGMAT)
        ALLOCATE(CGMAT(LMXX,LMXX,LMXX))
        CALL CLBSCH(LMXX,LMAX,CGMAT)
      END IF
      CG=CGMAT(LM1,LM2,LM3)
      RETURN
      END
!
!     ....................................................CLBSCH........
      SUBROUTINE CLBSCH(LMXX,LMAX,CG)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE CLEBSCH GORDAN COEFFICIENTS FOR REAL SPHERICAL    **
!     **  HARMONICS.                                                  **
!     **                                                              **
!     **    Y(LM1)*Y(LM2)= SUM(LM3|CG(LM1,LM2,LM3)*Y(LM3))            **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),INTENT(IN) :: LMXX
      INTEGER(4),INTENT(IN) :: LMAX
      REAL(8)   ,INTENT(OUT):: CG(LMXX,LMXX,LMXX)
      REAL(8)               :: FACT(0:50)
      REAL(8)               :: SQFACT(0:50)
      REAL(8)               :: PI,SQPI,SQ2
      INTEGER(4)            :: I,IMAX
      INTEGER(4)            :: LM1,L1,M1,N1,IS1
      INTEGER(4)            :: LM2,L2,M2,N2,IS2
      INTEGER(4)            :: LM3,L3,M3,N3,IS3
      INTEGER(4)            :: L3MIN,L3MAX
      REAL(8)               :: FAC0,FAC1,FAC2,SVAR
      INTEGER(4)            :: LMUP
      INTERFACE
        DOUBLE PRECISION FUNCTION THREEJ(L1,M1,L2,M2,L3,M3,FACT,SQFACT)
        INTEGER(4),INTENT(IN) :: L1,M1,L2,M2,L3,M3
        REAL(8)   ,INTENT(IN) :: FACT(0:50),SQFACT(0:50)
        END
      END INTERFACE
      INTERFACE
        DOUBLE PRECISION FUNCTION PRECG(L1,L2,L3,FACT,SQFACT)
        INTEGER(4),INTENT(IN) :: L1,L2,L3
        REAL(8)   ,INTENT(IN) :: FACT(0:50),SQFACT(0:50)
        END
      END INTERFACE
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      SQ2=DSQRT(2.D0)
      SQPI=DSQRT(PI)
      IMAX=4*LMAX+1
      IF(IMAX.GT.50) THEN
        CALL ERROR$MSG('DIMENSIONS IN CLBSCH TO SMALL')
        CALL ERROR$STOP('CLBSCH')
      END IF
      FACT(0)=1.D0
      SQFACT(0)=1.D0
      DO I=1,IMAX
        FACT(I)=FACT(I-1)*DBLE(I)
        SQFACT(I)=DSQRT(FACT(I))
      ENDDO
!
!     ==================================================================
!     ==  INITIALIZE ARRAY                                            ==
!     ==================================================================
      DO LM3=1,LMXX
        DO LM2=1,LMXX
          DO LM1=1,LMXX
            CG(LM1,LM2,LM3)=0.D0
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  CALCULATE CLEBSCH GORDAN COEFFICIENTS                       ==
!     ==================================================================
      DO L1=0,LMAX
        DO M1=-L1,L1
          N1=IABS(M1)
          IS1=0
          IF(M1.LT.0) IS1=1
!       ==
          DO L2=0,L1
            DO M2=-L2,L2
              N2=IABS(M2)
              IS2=0
              IF(M2.LT.0) IS2=1
!             ==
              L3MIN=IABS(L2-L1)
              L3MAX=MIN0(LMAX,L1+L2)
              DO L3=L3MIN,L3MAX,2
                DO M3=-L3,L3
                  N3=IABS(M3)
                  IS3=0
                  IF(M3.LT.0) IS3=1
!                 ==
                  IF(M1*M2.LT.0.AND.(M3.NE.-IABS(N1+N2) &
     &                         .AND. M3.NE.-IABS(N1-N2))) GOTO 210
                  IF(M1*M2.EQ.0.AND. M3.NE.M1+M2        ) GOTO 210
                  IF(M1*M2.GT.0.AND.(M3.NE.+IABS(N1+N2) &
     &                         .AND. M3.NE.+IABS(N1-N2))) GOTO 210
!                 ==
                  FAC0=PRECG(L1,L2,L3,FACT,SQFACT)
                  FAC1=0.D0
               IF(N1.EQ.0.AND.N2.EQ.0)FAC1=FAC1+FAC0/DSQRT(DBLE(2*L3+1))
                  FAC1=FAC1+THREEJ(L1, N1,L2, N2,L3, N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N3+IS3) &
     &                     +THREEJ(L1, N1,L2,-N2,L3,-N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N2+IS2) &
     &                     +THREEJ(L1,-N1,L2, N2,L3,-N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N1+IS1)
                  FAC2=DSQRT( 0.25D0*DBLE((2*L1+1)*(2*L2+1)) ) / SQ2 &
     &                          *(-1.D0)**(N3+(IS1+IS2+IS3)/2)    !WARNING CAN IS1+IS2+IS3 BE EVEN?
                  IF(M1.EQ.0)FAC2=FAC2/SQ2
                  IF(M2.EQ.0)FAC2=FAC2/SQ2
                  IF(M3.EQ.0)FAC2=FAC2/SQ2
                  LM1=(L1+1)**2-L1-M1
                  LM2=(L2+1)**2-L2-M2
                  LM3=(L3+1)**2-L3-M3
                  CG(LM1,LM2,LM3)=FAC0*FAC1*FAC2/SQPI
                  CG(LM2,LM1,LM3)=CG(LM1,LM2,LM3)
210               CONTINUE
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      IF(.NOT.TPR) RETURN
!     ==================================================================
!     == PRINTOUT FOR CHECK                                           ==
!     ==================================================================
      LMUP=(LMAX+1)**2
      DO LM1=1,LMUP
        DO LM2=1,LMUP
          DO LM3=1,LMUP
            IF(DABS(CG(LM1,LM2,LM3)).GT.1.D-6) THEN
              SVAR=CG(LM1,LM2,LM3)*DSQRT(4.D0*PI)
              WRITE(*,6000)LM1,LM2,LM3,CG(LM1,LM2,LM3),SVAR
6000          FORMAT(' CG ',3I3,2F10.5)
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................THREEJ..........
      DOUBLE PRECISION FUNCTION THREEJ(L1,M1,L2,M2,L3,M3,FACT,SQFACT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE WIEGNER 3J SYMBOLS                                **
!     **                                                              **
!     **  THREEJ = (-1)**L1-L2-M3 * ( L1  L2  L3 )                    **
!     **                            ( M1  M2  M3 )                    **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L1
      INTEGER(4),INTENT(IN) :: M1
      INTEGER(4),INTENT(IN) :: L2
      INTEGER(4),INTENT(IN) :: M2
      INTEGER(4),INTENT(IN) :: L3
      INTEGER(4),INTENT(IN) :: M3
      REAL(8)   ,INTENT(IN) :: FACT(0:50)
      REAL(8)   ,INTENT(IN) :: SQFACT(0:50)
      REAL(8)               :: R1,R2
      INTEGER(4)            :: KMAX,KMIN,K
      REAL(8)               :: SUM,SVAR
!     ******************************************************************
      IF(M3.NE.M1+M2) THEN
        THREEJ=0.D0
        RETURN
      END IF
!
!     ==================================================================
!     ==  CALCULATE THE 3J SYMBOLS  TIMES ....                        ==
!     ==================================================================
      R1=SQFACT(L1+L2-L3)*SQFACT(L3+L1-L2)*SQFACT(L3+L2-L1) &
     &  /SQFACT(L1+L2+L3+1)
      R2=SQFACT(L1+M1)*SQFACT(L1-M1)*SQFACT(L2+M2)*SQFACT(L2-M2) &
     &                              *SQFACT(L3+M3)*SQFACT(L3-M3)
      KMAX=MIN0(L1+L2-L3,L1-M1,L2+M2)
      KMIN=MAX0(0,L2-L3-M1,L1-L3+M2)
      SUM=0.D0
      DO K=KMIN,KMAX
        SVAR=FACT(K)*FACT(L1+L2-L3-K)*FACT(L1-M1-K)*FACT(L2+M2-K) &
     &      *FACT(L3-L2+M1+K)*FACT(L3-L1-M2+K)
        SUM=SUM+(-1.D0)**K/SVAR
      ENDDO
      THREEJ=R1*R2*SUM
      RETURN
      END
!
!     ..................................................PRECG...........
      DOUBLE PRECISION FUNCTION PRECG(L1,L2,L3,FACT,SQFACT)
!     ******************************************************************
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L1
      INTEGER(4),INTENT(IN) :: L2
      INTEGER(4),INTENT(IN) :: L3
      REAL(8)   ,INTENT(IN) :: FACT(0:50)
      REAL(8)   ,INTENT(IN) :: SQFACT(0:50)
      INTEGER(4)            :: LT,LTH
!     ******************************************************************
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      LT=L1+L2+L3
      LTH=LT/2
      IF(2*LTH.NE.LT) THEN
        PRECG=0.D0
      ELSE
        PRECG=DSQRT(DBLE(2*L3+1)/DBLE(LT+1))*FACT(LTH)/SQFACT(LT) &
     &       *SQFACT(LT-2*L1)/FACT(LTH-L1) &
     &       *SQFACT(LT-2*L2)/FACT(LTH-L2) &
     &       *SQFACT(LT-2*L3)/FACT(LTH-L3) &
     &       *(-1.D0)**(LTH-L3)
      END IF
      RETURN
      END
!
!.......................................................................
MODULE SPHERICALCCMAT_MODULE
!***********************************************************************
!**                                                                   **
!**  USED TO EVALUATE THE MATRIXES A0 AND A1 TO EVALUATE STRESSES     **
!**  IN A SPHERICAL HARMONICS REPRESENTATION                          **
!**                                                                   **
!**  GI(D/DGJ)SUM_L V_L Y_L                                           **
!**           = G_IG_J SUM_L [1/G^2 DV_L/DG -L V_L] Y_L               **
!**           + SUM_L SUM_L' V_L [P0_IJ(L,L')+PM_IJ(L,L')] Y_L'       **
!**                                                                   **
!**  WHERE P_0 IS NON-ZERO ONLY IF THE MAIN ANGULAR MOMENTA ARE       **
!**  IDENTICAL AND FOR P_M THE SECOND ANGULAR MOMENTUM IS LOWER       **
!**  BY 2 THAN THE FIRST                                              **
!**                                                                   **
!**  USE THE INTERFACES                                               **
!**    SPHERICAL$CCMAT0(LM1,I,LM2,CC)                                 **
!**    SPHERICAL$CCMATM(LM1,I,LM2,CC)                                 **
!**                                                                   **
!***********************************************************************
TYPE A_TYPE 
LOGICAL(4):: EXIST0(10)
LOGICAL(4):: EXISTM(10)
REAL(8)   :: A0(3,3,10)
REAL(8)   :: AM(3,3,10)
INTEGER(4):: LM20(10)
INTEGER(4):: LM2M(10)
END TYPE A_TYPE 
INTEGER(4)               :: LRXX=0
INTEGER(4)               :: LMRXX=0
TYPE(A_TYPE),ALLOCATABLE :: THIS(:)
END MODULE SPHERICALCCMAT_MODULE
!
!     ..................................................................
      SUBROUTINE SPHERICAL$CLEARCCMAT
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
!     *******************************************************************
      IF(ALLOCATED(THIS))DEALLOCATE(THIS)
      LRXX=0
      LMRXX=0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SPHERICALCCMAT_INITIALIZE(LM_)
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN):: LM_
      REAL(8)              :: A1(3,3)
      REAL(8)              :: FOURPIBY3
      REAL(8)              :: DSMALL=1.D-12
      REAL(8)              :: CG1,CG2
      INTEGER(4),PARAMETER :: MAPLM(3)=(/2,4,3/)
      INTEGER(4)           :: I,J,LMI,LMJ,LM1,LM2,LM3,LM3X,L1,L2,LM
      LOGICAL(4)           :: TCHK,T0,TM
!     ******************************************************************
      IF(ALLOCATED(THIS))DEALLOCATE(THIS)
      LRXX=MAX(INT(SQRT(REAL(LM_-1,KIND=8)+DSMALL)),2)
PRINT*,'INITIALIZE TO LRXX ',LRXX
      LMRXX=(LRXX+1)**2
      ALLOCATE(THIS(LMRXX))
      DO LM=1,LMRXX
        THIS(LM)%A0(:,:,:)=0.D0
        THIS(LM)%AM(:,:,:)=0.D0
        THIS(LM)%EXIST0(:)=.FALSE.
        THIS(LM)%EXISTM(:)=.FALSE.
        THIS(LM)%LM20(:)=0
        THIS(LM)%LM2M(:)=0
      ENDDO
      FOURPIBY3=16.D0*DATAN(1.D0)/3.D0  !4PI/3
!
!     ==================================================================
      DO LM1=1,LMRXX
        L1=INT(SQRT(REAL(LM1-1,KIND=8)+DSMALL))
        LM3X=L1**2
        DO LM2=1,LMRXX
          L2=INT(SQRT(REAL(LM2-1,KIND=8)+DSMALL))
          T0=L2.EQ.L1
          TM=L2.EQ.L1-2
          IF(.NOT.(T0.OR.TM)) CYCLE
!
!         ==============================================================
!         ==  EVALUATE MATRIX                                         ==
!         ==============================================================
          A1(:,:)=0.D0
          DO LM3=1,LM3X
            DO I=1,3
              LMI=MAPLM(I)
              CALL CLEBSCH(LMI,LM1,LM3,CG1)
              DO J=1,3
                LMJ=MAPLM(J)
                CALL CLEBSCH(LMJ,LM2,LM3,CG2)
                A1(I,J)=A1(I,J)+CG1*CG2
              ENDDO
            ENDDO
          ENDDO 
          A1(:,:)=FOURPIBY3*REAL(2*L1+1,KIND=8)*A1(:,:)
!
!         ==============================================================
!         ==  CYCLE IF THERE ARE ONLY ZERO ELEMENTS                    ==
!         ==============================================================
          TCHK=.FALSE.
          DO I=1,3
            DO J=1,3
              TCHK=TCHK.OR.(A1(I,J).NE.0.D0)
            ENDDO
          ENDDO
          IF(.NOT.TCHK) CYCLE                    
!
!         ==============================================================
!         ==  STORE L2=L1 ELEMENTS                                    ==
!         ==============================================================
          IF(T0) THEN
            I=1
            DO WHILE(THIS(LM1)%EXIST0(I)) 
              I=I+1
              IF(I.GT.10) THEN 
                CALL ERROR$MSG('SUMRULE VIOLATED - OVERFLOW A0')
                CALL ERROR$STOP('SPHERICAL_CCMAT')
              END IF
            END DO
            THIS(LM1)%A0(:,:,I)=A1(:,:)
            THIS(LM1)%LM20(I)=LM2 
            THIS(LM1)%EXIST0(I)=.TRUE.
          END IF
!
!         ==============================================================
!         ==  STORE L2=L1-2 ELEMENTS                                  ==
!         ==============================================================
          IF(TM) THEN
            I=1
            DO WHILE(THIS(LM1)%EXISTM(I)) 
              I=I+1
              IF(I.GT.10) THEN 
                CALL ERROR$MSG('SUMRULE VIOLATED - OVERFLOW AM')
                CALL ERROR$STOP('SPHERICAL_CCMAT')
              END IF
            END DO
            THIS(LM1)%AM(:,:,I)=A1(:,:)
            THIS(LM1)%LM2M(I)=LM2
            THIS(LM1)%EXISTM(I)=.TRUE.
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SPHERICAL$CCMAT0(LM1,I,LM2,CC)
!     ******************************************************************
!     **  MATRIX P0 USED TO EVALUATE THE STRESS IN A REPRESENTATION   **
!     **  OF SPHERICAL HARMONICS                                      **
!     **  SEE SPHERICALCCMAT_MODULE FOR DETAILS                       **
!     **  LM2=0 AND CC=0.D0 ON RETURN IF I IS TOO LARGE               **
!     ******************************************************************
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1     ! FIRST ANGULAR MOMENTUM INDEX  
      INTEGER(4),INTENT(IN) :: I       ! COUNTER FOR LM2               
      INTEGER(4),INTENT(OUT):: LM2     ! SECOND ANGULAR MOMENTUM INDEX 
      REAL(8)   ,INTENT(OUT):: CC(3,3) ! P0_IJ(LM1,LM2)                
      LOGICAL(4)            :: TCHK
!     ******************************************************************
      IF(LM1.GT.LMRXX) CALL SPHERICALCCMAT_INITIALIZE(LM1)
      TCHK=I.LE.10
      IF(TCHK) TCHK=THIS(LM1)%EXIST0(I) 
      IF(TCHK) THEN
        LM2=THIS(LM1)%LM20(I)
        CC=THIS(LM1)%A0(:,:,I)
      ELSE
        LM2=0
        CC(:,:)=0.D0
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SPHERICAL$CCMATM(LM1,I,LM2,CC)
!     ******************************************************************
!     **  MATRIX P0 USED TO EVALUATE THE STRESS IN A REPRESENTATION   **
!     **  OF SPHERICAL HARMONICS                                      **
!     **  SEE SPHERICALCCMAT_MODULE FOR DETAILS                       **
!     **  LM2=0 AND CC=0.D0 ON RETURN IF I IS TOO LARGE               **
!     ******************************************************************
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1      ! FIRST ANGULAR MOMENTUM INDEX  
      INTEGER(4),INTENT(IN) :: I        ! COUNTER FOR LM2               
      INTEGER(4),INTENT(OUT):: LM2      ! SECOND ANGULAR MOMENTUM INDEX 
      REAL(8)   ,INTENT(OUT):: CC(3,3)  ! PM_IJ(LM1,LM2)                
      LOGICAL(4)            :: TCHK
!     ******************************************************************
      IF(LM1.GT.LMRXX) CALL SPHERICALCCMAT_INITIALIZE(LM1)
      TCHK=I.LE.10
      IF(TCHK) TCHK=THIS(LM1)%EXISTM(I)
      IF(TCHK) THEN
        LM2=THIS(LM1)%LM2M(I)
        CC=THIS(LM1)%AM(:,:,I)
      ELSE
        LM2=0
        CC(:,:)=0.D0
      END IF
      RETURN
      END

