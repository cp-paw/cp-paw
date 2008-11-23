!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDBESSEL(R,K2,LMX,J)
!     **                                                                      **
!     ** constructs the regular solutions of the Helmholtz equation           **
!     **                   (nabla^2 + k2)*psi(r)=0                            **
!     **                                                                      **
!     ** The solution behaves at the origin like                              **
!     **    J_{l,m}(r)=1/(2*l+1)!! * |r|^l  Y_{l,m}(r)  +O(r^l+1)             **
!     **                                                                      **
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(out):: J(LMX)
      INTEGER(4)            :: LX
      REAL(8)               :: K 
      REAL(8)               :: X,Y,dydx
      INTEGER(4)            :: LM,m
      INTEGER(4)            :: L
      REAL(8)               :: PI
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      CALL SPHERICAL$YLM(LMX,R,J)
      LX=int(SQRT(REAL(LMX+1.D-5)))-1
      K=SQRT(ABS(K2))     
      X=SQRT(SUM(R**2))
      LM=0
      DO L=0,LX
        IF(K2.GT.0.D0) THEN
          CALL SPFUNCTION$BESSEL(L,k*X,Y,dydx)  ! ABRAMOWITZ 10.1.25
          Y=Y/K**L
        ELSE IF(K2.LT.0.D0) THEN
          CALL SPFUNCTION$MODBESSEL(L,k*X,Y,dydx) !ABRAMOWITZ 10.2.4
          Y=Y/K**L
        ELSE
          CALL SPFUNCTION$BESSEL0(L,X,Y,dydx)  ! ABRAMOWITZ 10.1.2
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
      SUBROUTINE LMTO$SOLIDHANKEL(R,rad,K2,LMX,H)
!     **                                                                      **
!     ** constructs the irregular solutions of the Helmholtz equation         **
!     **                   (nabla^2 + k2)*psi(r)=0                            **
!     **                                                                      **
!     ** The solution behaves at the origin like                              **
!     **    H_{l,m}(r)=1/(2*l-1)!! * |r|^{-l-1}  Y_{l,m}(r)                   **
!     **                                                                      **
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R(3) ! position relative to the origin
      real(8)   ,intent(in) :: rad  ! inside rad the function is regularized
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(out):: H(LMX)  ! solid hankel function
      INTEGER(4)            :: LX
      REAL(8)               :: K 
      REAL(8)               :: X,xr,Y,dydx
      INTEGER(4)            :: LM,m
      INTEGER(4)            :: L
      REAL(8)               :: PI
      logical(4)            :: tcap
      real(8)               :: a,b
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      CALL SPHERICAL$YLM(LMX,R,H)
      LX=int(SQRT(REAL(LMX+1.D-5)))-1
      K=SQRT(ABS(K2))     
      Xr=SQRT(SUM(R**2))
      tcap=xr.lt.rad
      x=max(xr,rad)
      LM=0
      DO L=0,LX
        IF(K2.GT.0.D0) THEN
          CALL SPFUNCTION$NEUMANN(L,k*X,Y,dydx)  ! ABRAMOWITZ 10.1.26
          Y=-Y*K**(L+1)
          dYdx=-dYdx*K**(L+2)
        ELSE IF(K2.LT.0.D0) THEN
          CALL SPFUNCTION$MODHANKEL(L,k*X,Y,dydx) !ABRAMOWITZ 10.2.4
          Y=Y*2.D0/PI*K**(L+1)
          dYdx=dYdx*2.D0/PI*K**(L+2)
        ELSE
!         ==  y(x)= 1/(2l-1)!! * x**(-l-1) 
          CALL SPFUNCTION$NEUMANN0(L,X,Y,dydx)  ! ABRAMOWITZ 10.2.5
          Y=-Y     !
          dydx=-dydx 
        END IF 
!
!       == inside rad, match a parabola times r**l ==============================
        if(tcap) then
          b=0.5d0*(dydx*x-real(l,kind=8)*y)/x**(l+2)
          a=y/x**l-b*x**2
          if(l.eq.0) then
            y=a+b*xr**2
            dydx=2.d0*b*xr
          else
            y=a*xr**l+b*xr**(l+2)
            dydx=real(l,kind=8)*a*xr**(l-1)+real(l+2)*b*xr**(l+1)
          end if
        end if  
        DO M=1,2*L+1
          LM=LM+1
          H(LM)=H(LM)*Y
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      subroutine lmto$kbarmulticenter(n,norb,qbar,sbar,c)
!     **************************************************************************
!     ** determines the coefficients for the multicenter expansion of the     **
!     ** screened Hankelfunctions in termes of unscreened ones.               **
!     **************************************************************************
      implicit none
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(IN) :: QBAR(n)
      REAL(8)   ,INTENT(in) :: SBAR(N,NORB)
      REAL(8)   ,INTENT(out):: c(N,NORB)
      integer(4)            :: i
!     **************************************************************************
!
!     ==========================================================================
!     == SET UP coefficients for the multicenter expansion of kbar            ==
!     ==========================================================================
      do i=1,norb
        c(:,i)=qbar(:)*sbar(:,i)
        c(i,i)=c(i,i)+1.d0
      enddo
      return
      end subroutine lmto$kbarmulticenter

!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$expandqbar(NAT,lxx,lx,qbar,n,qbarvec)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: Nat
      INTEGER(4),INTENT(IN) :: LXx
      INTEGER(4),INTENT(IN) :: LX(NAT)
      REAL(8)   ,INTENT(IN) :: qbar(lxx+1,Nat)
      INTEGER(4),INTENT(IN) :: n
      REAL(8)   ,INTENT(out):: qbarvec(n)
      integer(4)            :: i,iat,l,im
!     **************************************************************************
      if(n.ne.sum((lx+1)**2)) then
        call error$stop('lmto$expandqbar')
      end if
      i=0
      do iat=1,nat
        do l=0,lx(iat)
          do im=1,2*l+1
            i=i+1
            qbarvec(i)=qbar(l+1,iat)
          enddo
        enddo
      enddo
      return
      end

!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SCREENEDSOLIDHANKEL(K2,NAT,RPOS,RAD,LX,N,NORB,C,R,F)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: K2
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RPOS(3,NAT)
      REAL(8)   ,INTENT(IN) :: Rad(NAT)
      INTEGER(4),INTENT(IN) :: LX(NAT)
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: Norb
      REAL(8)   ,INTENT(IN) :: C(N,norb)
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(out):: F(norb)
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
        CALL LMTO$SOLIDHANKEL(R(:)-RPOS(:,IAT),rad(iat),K2,LMX,H)
        F(:)=F(:)+matmul(H(:LMX),C(I+1:I+LMX,:))
        I=I+LMX
      ENDDO
      DEALLOCATE(H)
      RETURN
      END SUBROUTINE LMTO$SCREENEDSOLIDHANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDBESSELrad(l,R,K2,J,jder)
!     **                                                                      **
!     ** constructs the radial part of the regular solutions                  **
!     **  of the Helmholtz equation  (nabla^2 + k2)*psi(r)=0                  **
!     **                                                                      **
!     ** The solution behaves at the origin like                              **
!     **    J_{l,m}(r)=1/(2*l+1)!! * |r|^l  Y_{l,m}(r)  +O(r^l+1)             **
!     **                                                                      **
      IMPLICIT NONE
      integer(4),intent(in) :: l
      REAL(8)   ,INTENT(IN) :: R
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      REAL(8)   ,INTENT(out):: J
      REAL(8)   ,INTENT(out):: Jder
      REAL(8)               :: K 
      REAL(8)               :: X,Y,dydx
      REAL(8)               :: PI
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      K=SQRT(ABS(K2))     
      X=r
      IF(K2.GT.0.D0) THEN
        CALL SPFUNCTION$BESSEL(L,k*X,Y,dydx)  ! ABRAMOWITZ 10.1.25
        Y=Y/K**L
        dydx=dydx/K**(L-1)
      ELSE IF(K2.LT.0.D0) THEN
        CALL SPFUNCTION$MODBESSEL(L,k*X,Y,dydx) !ABRAMOWITZ 10.2.4
        Y=Y/K**L
        dYdx=dYdx/K**(L-1)
      ELSE
        CALL SPFUNCTION$BESSEL0(L,X,Y,dydx)  ! ABRAMOWITZ 10.1.2
      END IF
      j=y
      jder=dydx
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDHANKELrad(l,R,K2,H,hder)
!     **                                                                      **
!     ** constructs radial part of the the irregular solutions                **
!     ** of the Helmholtz equation (nabla^2 + k2)*psi(r)=0                    **
!     **                                                                      **
!     ** The solution behaves at the origin like                              **
!     **    H_{l,m}(r)=1/(2*l-1)!! * |r|^{-l-1}  Y_{l,m}(r)                   **
!     **                                                                      **
      IMPLICIT NONE
      integer(4),INTENT(IN) :: l    ! main angular momentum
      REAL(8)   ,INTENT(IN) :: R    ! radius
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      REAL(8)   ,INTENT(out):: H    ! radial part of the hankel function
      REAL(8)   ,INTENT(out):: Hder ! radial derivative of the hankel function
      REAL(8)               :: K 
      REAL(8)               :: X,Y,dydx
      REAL(8)               :: PI
      REAL(8)               :: svar
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      K=SQRT(ABS(K2))     
      X=r
      IF(K2.GT.0.D0) THEN
        CALL SPFUNCTION$NEUMANN(L,k*X,Y,dydx)  ! ABRAMOWITZ 10.1.26
        svar=-K**(L+1)
        Y=svar*Y
        dydx=svar*dydx*k
      ELSE IF(K2.LT.0.D0) THEN
        CALL SPFUNCTION$MODHANKEL(L,k*X,Y,dydx) !ABRAMOWITZ 10.2.4
        svar=2.D0/PI*K**(L+1)
        Y=svar*y
        dydx=svar*dydx*k
      ELSE
!       ==  y(x)= 1/(2l-1)!! * x**(-l-1) 
        CALL SPFUNCTION$NEUMANN0(L,X,Y,dydx)  ! ABRAMOWITZ 10.2.5
        Y=-Y     ! 
        dydx=-dydx
      END IF   
      h=y
      hder=dydx
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$STRUCTURECONSTANTS(R1,L1X,R2,L2X,K2,S)
!     **                                                                      **
!     **  constructs the structure constants that mediate an expansion        **
!     **  of a solid hankel function H_{l,m}(r-R1) centered at R1             **
!     **  into solid bessel functions  J_{l,m}(r-r2) centered at r2           **
!     **                                                                      **
!     **    H_{l,m}(r-r1) = sum_{l',m'} J_{l',m'}(r-r2) * S_{l',m',l,m}       **
!     **                                                                      **
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1(3)
      REAL(8)   ,INTENT(IN) :: R2(3)
      INTEGER(4),INTENT(IN) :: L1X
      INTEGER(4),INTENT(IN) :: L2X
      REAL(8)   ,INTENT(IN) :: K2 ! 2ME/HBAR**2
      REAL(8)   ,INTENT(OUT):: S((L2X+1)**2,(L1X+1)**2)
      real(8)   ,parameter  :: rad=1.d-6
      REAL(8)               :: k
      REAL(8)               :: PI
      REAL(8)               :: SVAR
      INTEGER(4)            :: L3X,LM1X,LM2X,LM3X
      INTEGER(4)            :: LM1,LM2,LM3,l,l1,l2,l3,IM,LM3A,LM3B
      INTEGER(4)            :: LOFLM((L1X+L2X+1)**2)
      REAL(8)               :: H((L1X+L2X+1)**2)
      REAL(8)               :: HANKEL ! HANKEL FUNCTION OF THE DISTANCE
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      complex(8)            :: kappa
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      L3X=L1X+L2X
      LM1X=(L1X+1)**2
      LM2X=(L2X+1)**2
      LM3X=(L3X+1)**2
      lm3=0
      do l=0,l3x
        do Im=1,2*l+1
          lm3=lm3+1
          loflm(lm3)=L
        ENDDO
      ENDDO
      if(k2.gt.0.d0) then
        kappa=cmplx(0.d0,-sqrt(k2))
      else if(k2.lt.0.d0) then
        kappa=cmplx(sqrt(-k2),0.d0)
      else
        kappa=(0.d0,0.d0)
      end if

!     == CALCULATE HANKEL FUNCTION OF THE DISTANCE =============================
      call LMTO$SOLIDHANKEL(R2-r1,rad,K2,LM3X,H)
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
          DO LM3=1,lm3x
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
!     == MULTIPLY WITH -4*PI (-1)**l2 ============================================
      S=-4.D0*PI*S
      do lm2=1,lm2x
        s(lm2,:)=s(lm2,:)*(-1.d0)**loflm(lm2)
      enddo
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$Q(L,RAD,VAL,DER,K2,QPAR)
!     **************************************************************************      
!     ** parameter needed to  screen the structure constants                  **
!     **   |K>-|J>qbar has the same logarithmic derivative as |phidot>.       **
!     ** val and der are value and derivative of phidot.                       **
!     **                                                                      **
!     **************************************************************************      
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: VAL
      REAL(8)   ,INTENT(IN) :: DER
      REAL(8)   ,INTENT(IN) :: K2
      REAL(8)   ,INTENT(OUT):: QPAR
      REAL(8)               :: Jval,JDER
      REAL(8)               :: Kval,KDER
!     **************************************************************************      
      CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,Jval,JDER)
      CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,Kval,KDER)
      QPAR=(JVAL*DER-VAL*JDER)/(KVAL*DER-VAL*KDER)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SCREEN(TSTART,N,NORB,QBAR,S,SBAR)
!     **************************************************************************      
!     **  determines screened structure constants for a cluster               **
!     **      |Kbar_i>=sum_j |K_j> (delta_ji+Qbar_j*sbar_ji)                  **
!     **  start with sbar=0 or give better estimate                           **
!     **                                                                      **
!     **                                                                      **
!     **      sbar=                                                           **
!     **                                                                      **
!     **************************************************************************      
      IMPLICIT NONE
      logical(4),intent(in) :: tstart
      INTEGER(4),INTENT(IN) :: N            ! #(orbitals on the cluster)
      INTEGER(4),INTENT(IN) :: NORB         ! #(orbitals on central site
      REAL(8)   ,INTENT(IN) :: QBAR(N)      !
      REAL(8)   ,INTENT(IN) :: S(N,N)       ! unscreened-structure constants
      REAL(8)   ,INTENT(inOUT):: SBAR(N,NORB) !screened structure constants
      real(8)   ,parameter  :: tol=1.d-5    ! tolerance for converence
      integer(4),parameter  :: niter=1000   ! x#(iterations)
      real(8)               :: alpha=0.5d0  ! mixing factor
      real(8)               :: dsbar(n,norb)
      real(8)               :: s0(n,norb)
      real(8)               :: a(n,n)
      integer(4)            :: i
      real(8)               :: delta
      integer(4)            :: iter
      logical(4)            :: convg
!     **************************************************************************      
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      if(tstart) then
        do i=1,n
          a(:,i)=-qbar(:)*s(:,i)
          a(i,i)=a(i,i)+1.d0
        enddo
        s0(:,:)=s(:,:norb)
        call LIB$MATRIXSOLVER8(N,n,norb,A,sbar,s0)
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      else
        s0(:,:)=s(:,:norb)
        do iter=1,niter
          do i=1,norb
            dsbar(:,i)=qbar(:)*sbar(:,i)
            dsbar(i,i)=1.d0+dsbar(i,i)
          enddo
          dsbar=matmul(s,dsbar)-sbar
          delta=maxval(abs(dsbar))
print*,'lmto$screen: delta ',iter,delta
          convg=delta.lt.tol
          if(convg) exit
          sbar=sbar+dsbar*alpha
        enddo
        if(.not.convg) then
          call error$msg('loop not converged')
          call error$stop('lmto$screen')
        end if
      end if
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT,RPOS,LX,QbAR,N,NORB,SBAR)
!     **************************************************************************      
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **                                                                      **
!     ** remark: ThE CENTRAL ATOM IS THE FIRST ATOM IN THE LIST               **
!     **                                                                      **
!     ** remark: Initialize sbar with zero or a better estimate               **
!     **                                                                      **
!     **                                                                      **
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **************************************************************************      
      implicit none
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RPOS(3,NAT)
      INTEGER(4),INTENT(IN) :: LX(NAT)
      REAL(8)   ,INTENT(IN) :: QBAR(n)
      REAL(8)   ,INTENT(IN) :: K2
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(inOUT):: SBAR(N,NORB)
      INTEGER(4)            :: II,IAT,LN,L,IM
      INTEGER(4)            :: I1,I2
      INTEGER(4)            :: iat1,iat2
      REAL(8)               :: R1(3),R2(3)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: L1X,L2X
      REAL(8)               :: S0(n,n)
      REAL(8)  ,allocatable :: S1(:,:)
      REAL(8)               :: qbarvec(n)
!     **************************************************************************      
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
          if(iat2.eq.iat1) then
            I2=I2+lmn1
            cycle
          end if
          R2(:)=RPOS(:,IAT2)
write(*,fmt='("r1=",3f10.4," dr=",3f10.4)')r1,r2-r1
          L2X=LX(IAT2)
          LMN2=(L2X+1)**2
          CALL LMTO$STRUCTURECONSTANTS(R1,L1X,R2,L2X,K2,S1(:lmn2,:))
          S0(I2+1:I2+LMN2,i1+1:I1+LMN1)=S1(:lmn2,:)
          I2=I2+(L2X+1)**2
        ENDDO
        DEALLOCATE(S1)
        I1=I1+(L1X+1)**2
      ENDDO
!
!     ==========================================================================
!     == SCREEN STRUCTURE CONSTANTS                                           ==
!     ==========================================================================
      sbar=0.d0
      CALL LMTO$SCREEN(.false.,N,NORB,QBAR,S0,SBAR)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$neighborlist(rbas,nat,r,rc,nnx,NNB,nnlist)
!     **************************************************************************
!     **  this is a simple neighborlist routine                               **
!     **************************************************************************
      implicit none
      integer(4)   ,intent(in) :: nat       ! #(atoms)
      real(8)      ,intent(in) :: rbas(3,3) ! lattice vectors
      real(8)      ,intent(in) :: r(3,nat)  ! atom positions
      real(8)      ,intent(in) :: rc        ! cutoff radius for the neighorlist
      integer(4)   ,intent(in) :: nnx       ! x#(neighbors per atom)
      INTEGER(4)   ,intent(out):: nnB       ! #(NEIGHBORS)
      INTEGER(4)   ,intent(out):: nnlist(5,nnx) ! neighborlist (IAT1,IAT2,IT(3))
      real(8)                  :: rbasinv(3,3)
      real(8)                  :: rfold(3,nat)
      real(8)                  :: x(3)      ! relative coordinates
      real(8)                  :: xfold(3)      ! relative coordinates
      integer(4)               :: itfold(3,nat)   ! shift
      integer(4)               :: iat,i,iat1,iat2,it1,it2,it3
      real(8)                  :: rmax2     ! squared cutoff radius
      real(8)                  :: tvec(3)
      integer(4)               :: itvec(3)
      integer(4)               :: min1,max1,min2,max2,min3,max3
      real(8)                  :: d(3),d2
      real(8)                  :: x0,y0,z0
!     **************************************************************************
!
!     ==========================================================================
!     == fold atom positions into the first unit cell                         ==
!     ==========================================================================
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      DO IAT=1,NAT
        X(:)=MATMUL(RBASINV,R(:,IAT))
        DO I=1,3
          Xfold(I)=MODULO(X(I),1.D0)
        ENDDO
        itfold(:,iat)=nint(xfold-x)
        RFOLD(:,IAT)=MATMUL(RBAS,Xfold)
      ENDDO
!
!     ==========================================================================
!     == fold atom positions into the first unit cell                         ==
!     ==========================================================================
      rmax2=rc**2
      nnb=0
      DO IAT1=1,NAT
!       == place onsite element for each atom first in the neighborlist
        nnb=nnb+1
        IF(NNB.gt.NNX) THEN
          call error$msg('maximum number of neighbors exceeded')
          call error$i4val('nnb',nnb)
          call error$i4val('nnx',nnx)
          call error$stop('lmto$neighborlist')
        end if
        nnlist(1,nnb)=iat1
        nnlist(2,nnb)=iat1
        nnlist(3:5,nnb)=0
        do iat2=1,nat
          X0=Rfold(1,IAT1)-rfold(1,iat2)
          Y0=Rfold(2,IAT1)-rfold(2,iat2)
          Z0=Rfold(3,IAT1)-rfold(3,iat2)
          CALL BOXSPH(rbas,X0,Y0,Z0,RC,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!         == LOOP OVER BOXES IN THE NEIGHBORHOOD ================================
          DO It1=MIN1,MAX1
            DO It2=MIN2,MAX2
              DO It3=MIN3,MAX3
                if(iat1.eq.iat2.and.it1.eq.0.and.it2.eq.0.and.it3.eq.0) cycle
                itvec(1)=it1
                itvec(2)=it2
                itvec(3)=it3
                tvec(:)=matmul(rbas,real(itvec,kind=8))
!               == distance criterion =========================================
                D(:)=Rfold(:,IAT2)+tvec(:)-Rfold(:,IAT1)
                D2=SUM(D(:)**2)
                IF(D2.GT.RMAX2) CYCLE
                NNB=NNB+1
                IF(NNB.gt.NNX) THEN
                  call error$msg('maximum number of neighbors exceeded')
                  call error$i4val('nnb',nnb)
                  call error$i4val('nnx',nnx)
                  call error$stop('lmto$neighborlist')
                end if
                NnLIST(1,NNB)=IAT1
                NnLIST(2,NNB)=IAT2
                itvec(:)=itvec(:)+itfold(:,iat2)-itfold(:,iat1)
                NnLIST(3:5,NNB)=Itvec(:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!!$do iat1=1,nat
!!$ d(:)=r(:,iat1)-rfold(:,iat1)
!!$ write(*,fmt='("d ",3f10.5,3i5)')d(:),itfold(:,iat1)
!!$enddo
!!$do i=1,nnb
!!$  iat1=nnlist(1,i)
!!$  iat2=nnlist(2,i)
!!$  itvec(:)=nnlist(3:5,i)
!!$  d(:)=r(:,iat2)-r(:,iat1)+matmul(rbas,real(itvec,kind=8))
!!$  write(*,fmt='(i5," dis ",f10.5," d ",3f10.5)')i,sqrt(sum(d**2)),d(:)
!!$enddo
!!$print*,'rc ',rc
!!$stop
      return
      end
