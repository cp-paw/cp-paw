!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine specialfunction$test()
!      *************************************************************************
!      **                                                                     **
!      **                                                                     **
!      **                                                                     **
!      *************************************************************************
       implicit none
       integer(4),parameter :: lx=3
       integer(4),parameter :: nx=200
       integer(4),parameter :: nfil=101
       real(8)   ,parameter :: xmax=10.d0
       real(8)   ,parameter :: dx=xmax/(nx-1) ! rounding errors acceptable
       real(8)              :: x
       real(8)              :: y(lx+1)
       real(8)              :: dydx(lx+1)
       integer(4)           :: ix,l
!      *************************************************************************
!
!      =========================================================================
!      ==  bessel functions                                                   ==
!      =========================================================================
       open(nfil,file='test_bessel.dat')
       do ix=1,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$BESSEL(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  modified bessel functions                                          ==
!      =========================================================================
       open(nfil,file='test_modbessel.dat')
       do ix=1,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$modBESSEL(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  bessel functions for kappa=0                                       ==
!      =========================================================================
       open(nfil,file='test_bessel0.dat')
       do ix=1,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$BESSEL0(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  modified hankel functions                                          ==
!      =========================================================================
       open(nfil,file='test_modhankel.dat')
       do ix=2,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$modhankel(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  neumann functions                                                  ==
!      =========================================================================
       open(nfil,file='test_neumann.dat')
       do ix=2,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$neumann(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  modified neumann functions                                         ==
!      =========================================================================
       open(nfil,file='test_modneumann.dat')
       do ix=2,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$modneumann(L,X,Y(l+1),dydx(l+1))       
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
!
!      =========================================================================
!      ==  neumann functions for kappa=0                                      ==
!      =========================================================================
       open(nfil,file='test_neumann0.dat')
       do ix=2,nx
         x=dx*real(ix-1,kind=8)
         do l=0,lx
           call Spfunction$NEUMANN0(L,X,Y(l+1),dydx(l+1))
         enddo
         write(nfil,*)x,y,dydx
       enddo
       close(nfil)
       return
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine specialfunction$erf(x,val)
!      *************************************************************************
!      **  returns the error function erf(x)
!      ** see numerical recipes
!      *************************************************************************
       implicit none
       real(8),intent(in) :: x
       real(8),intent(out):: val
!      *************************************************************************
       call specialfunction$gammp(x**2,0.5d0,val)
       If(x.lt.0.d0) val=-val
       return
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine specialfunction$gammp(x,a,val)
!      *************************************************************************
!      **  returns the incomplete gamma function P(a,x)
!      **  see numerical recipes
!      *************************************************************************
       implicit none
       real(8),intent(in) :: x
       real(8),intent(in) :: a
       real(8),intent(out):: val
!      *************************************************************************
       if(x.lt.0.d0.or.a.le.0.d0) then
         stop 'error stop in specialfunction$gammp'
       end if
       if(x.lt.a+1.d0) then
!        == use series expansion ===================================
         call specialfunction_gser(x,a,val)
       else
!        == use continued fraction representation ==================
         call specialfunction_gcf(x,a,val)
         val=1.d0-val
       end if
       return
       end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_GAMMLN(X,VAL)
!      ** ln(gamma(xx))
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: X
       REAL(8),INTENT(OUT):: VAL
       REAL(8),PARAMETER  :: COF(6)=(/76.18009173D0,-86.50532033D0 &
      &                             ,24.01409822D0,-1.231739516D0 &
      &                             ,0.120858003D-2,-0.536382D-5/)
       REAL(8),PARAMETER  :: STP=2.50662827465D0
       real(8)            :: x1,tmp,ser
       integer(4)         :: j
!      *************************************************************
       X1=X-1.D0
       TMP=X1+5.5D0
       TMP=(X1+0.5D0)*LOG(TMP)-TMP
       SER=1.D0
       DO J=1,6
         X1=X1+1.D0
         SER=SER+COF(J)/X1
       ENDDO
       val=TMP+LOG(STP*SER)
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_Gser(X,A,val)
!      ** 
!      ** incomplete gamma function P(a,x) evaluated by its series expansion
!      ** 
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: X
       REAL(8)   ,INTENT(IN) :: A
       REAL(8)   ,INTENT(OUT):: VAL
       integer(4),parameter  :: itmax=100
       real(8)   ,parameter  :: eps=3.d-7
       real(8)               :: ap,sum,del,gln
       integer(4)            :: n
!      *************************************************************
       if(x.lt.0.d0) then
         stop 'invalid argument for gser'
       end if
       if(x.eq.0.d0) then
         val=0.d0
         return
       end if
       ap=a
       sum=1.d0/a
       del=sum
       do n=1,itmax
         ap=ap+1.d0
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*eps) then
           call specialfunction_gammln(a,gln)
           val=sum*exp(-x+a*log(x)-gln)
           return
         end if
       enddo
       stop 'a too large, itmax too small; stop in gser'
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_Gcf(X,A,val)
!      ** 
!      ** incomplete gamma function q(a,x) evaluated by its continued 
!      **  fraction representation
!      ** 
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: X
       REAL(8)   ,INTENT(IN) :: A
       REAL(8)   ,INTENT(OUT):: VAL
       integer(4),parameter  :: itmax=100
       real(8)   ,parameter  :: eps=3.d-7
       real(8)               :: gold,gln
       real(8)               :: a0,a1,b0,b1
       real(8)               :: fac,an,ana,anf,g
       integer(4)            :: n
!      *************************************************************
       call specialfunction_gammln(a,gln)
       gold=0.d0
       a0=1.d0
       a1=x
       b0=0.d0
       b1=1.d0
       fac=1.d0
       do n=1,itmax
         an=real(n,kind=8)
         ana=an-a
         a0=(a1+a0*ana)*fac
         b0=(b1+b0*ana)*fac
         anf=an*fac
         a1=x*a0+anf*a1
         b1=x*b0+anf*b1
         if(a1.ne.0.d0) then
           fac=1.d0/a1
           g=b1*fac
           if(abs((g-gold)/g).lt.eps) then
             val=exp(-x+a*log(x)-gln)*g
             return
            end if
            gold=g
          endif
        enddo
        stop 'a too large, itmax soo small, stop in gcf'
        end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPECIALFUNCTION$BESSEL(L,X,Y)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.2              FOR   X < L                           **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > L                           **
!     **                                                                      **
!     ** todo: this routine is duplicated by spfunction$bessel                **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)               :: TRIG(4)
      REAL(8)               :: FACUL(0:100)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: PI
      REAL(8)               :: ARG
      REAL(8)               :: XSQ
      REAL(8)               :: FAC
      INTEGER(4)            :: I,K,IL,II,ISVAR
!     ******************************************************************
      IF(X.GT.real(L,kind=8)) THEN
        PI=4.D0*ATAN(1.D0)
        ARG=X-0.5D0*real(L,kind=8)*PI
        TRIG(1)=SIN(ARG)/X
        TRIG(2)=COS(ARG)/X
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        Y=TRIG(1)
        IF(L.EQ.0) RETURN
!       ==  double faculty facul(l)=(2*l)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*real(I,kind=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
!         FAC=FAC*XSQ*DBLE(L+K)/DBLE(K*(L-K))
          Y=Y+FAC*TRIG(II)
        ENDDO
!       II=MOD(L,4)+1
!       FAC=FAC*XSQ*DBLE(2*L)/DBLE(L)
!       Y=Y+FAC*TRIG(II)
        RETURN
      END IF
!     ==================================================================
!     ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                        ==
!     ==================================================================
      ISVAR=1
      DO IL=1,L
        ISVAR=ISVAR*(2*IL+1)
      ENDDO
      IF(L.NE.0.D0) THEN
        FAC=X**L/DBLE(ISVAR)
      ELSE
        FAC=1.D0/DBLE(ISVAR)
      END IF
      Y=FAC
      XSQ=-0.5D0*X*X
      ISVAR=2*L+1
      DO I=1,1000
        ISVAR=ISVAR+2
        FAC=FAC*XSQ/DBLE(I*ISVAR)
        Y=Y+FAC
        IF(ABS(FAC).LT.TOL) GOTO 9999
      ENDDO
      CALL ERROR$MSG('Y NOT CONVERGED')
      CALL ERROR$I4VAL('L',L)
      CALL ERROR$R8VAL('X',X)
      CALL ERROR$STOP('specialfunction$bessel')
9999  CONTINUE
      RETURN
      END
!
!===============================================================================
!===============================================================================
!===============================================================================
!=====     The following routines need to be tested.                          ==
!===============================================================================
!===============================================================================
!===============================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Spfunction$BESSEL(L,X,Y,dydx)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.2              FOR   X < L                           **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > L                           **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 2009 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dYdx ! derivative
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: TRIG(4)
      REAL(8)               :: dTRIG(4)
      REAL(8)               :: FACUL(0:2*l)
      REAL(8)               :: PI
      REAL(8)               :: ARG
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,dfac
      INTEGER(4)            :: I,K,IL,II,ISVAR
      logical(4)            :: convg
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('bessel function FOR NEGATIVE ARGuMENTS undefined')
        CALL ERROR$STOP('SPFUNCTION$bessel')
      else IF(X.le.real(L,kind=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS  x<L                         ==
!       ========================================================================
        ISVAR=1
        DO IL=1,L
          ISVAR=ISVAR*(2*IL+1)
        ENDDO
        if(l.eq.0) then
          FAC=1.D0/real(ISVAR,kind=8)
          dfac=0.d0
        else if(l.eq.1) then
          FAC=X/real(ISVAR,kind=8)
          dfac=1.d0/real(ISVAR,kind=8)
        else 
          FAC=X**L/real(ISVAR,kind=8)
          dFAC=real(l,kind=8)*X**(L-1)/real(ISVAR,kind=8)
        END IF
        Y=FAC
        dYdx=dFAC
        XSQ=-0.5D0*X*X
        ISVAR=2*L+1
        DO I=1,1000
          ISVAR=ISVAR+2
!         = do derivative before value because value is modified ===============
          dFAC=(dFAC*XSQ-FAC*X)/real(I*ISVAR,kind=8)
          FAC=FAC*XSQ/real(I*ISVAR,kind=8)
          Y=Y+FAC
          dYdx=dYdx+dFAC
          convg=abs(FAC).LT.TOL
          if(convg) exit
        ENDDO
        if(.not.convg) then 
          CALL ERROR$MSG('Y NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('X',X)
          CALL ERROR$STOP('spfunction$bessel')
        end if
      else 
!       ========================================================================
!       ==  EXPANSION FOR large ARGUMENTS  x>L                                ==
!       ========================================================================
        PI=4.D0*ATAN(1.D0)
        ARG=X-0.5D0*real(L,kind=8)*PI
        TRIG(1)=SIN(ARG)/X
        TRIG(2)=COS(ARG)/X
        dtrig(1)=cos(arg)/x-sin(arg)/x**2
        dTRIG(2)=-sin(ARG)/X-COS(ARG)/X**2
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        DTRIG(3)=-DTRIG(1)
        DTRIG(4)=-DTRIG(2)
        Y=TRIG(1)
        DYDX=DTRIG(1)
        IF(L.EQ.0) RETURN
!       ==  DOUBLE FACULTY FACUL(L)=(2*L)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*REAL(I,KIND=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
          DFAC=-REAL(K,KIND=8)*FAC/X
!         FAC=FAC*XSQ*DBLE(L+K)/DBLE(K*(L-K))
          Y=Y+FAC*TRIG(II)
          DYDX=DYDX+DFAC*TRIG(II)+FAC*DTRIG(II)
        ENDDO
!       II=MOD(L,4)+1
!       FAC=FAC*XSQ*DBLE(2*L)/DBLE(L)
!       Y=Y+FAC*TRIG(II)
        RETURN
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$NEUMANN(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  SPHERICAL BESSEL FUNCTION OF THE SECOND KIND Y_L(X)                 **
!     **    Y = -X^L*(-1/X D/DX)^L [COS(X)/X]  FORMULA 10.1.26                **
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION                           **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.3   FOR   X < L                                      **
!     **    FORMULA 10.1.9   FOR   X > L                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE OF NEUMANN FUNCTION AT X
      REAL(8)               :: TRIG(4)
      REAL(8)               :: DTRIG(4)
      REAL(8)               :: FACUL(0:2*l)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: PI
      REAL(8)               :: ARG
      REAL(8)               :: XSQ,dxsq
      REAL(8)               :: FAC,dfac
      INTEGER(4)            :: I,K,IL,II,ISVAR
      real(8)               :: m1powerl
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION NOT DEEFINED FOR NEGATIVE ARGYMENTS')
        CALL ERROR$i4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.EQ.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION WITH ZERO ARGUMENT DIVERGES')
        CALL ERROR$i4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.LE.REAL(L,KIND=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS   x<l                        ==
!       ========================================================================
        ISVAR=-1
        DO IL=2,L
          ISVAR=ISVAR*(2*IL-1)
        ENDDO
        FAC=REAL(ISVAR,KIND=8)/X**(L+1)
        DFAC=-REAL(L+1,KIND=8)*FAC/X
        Y=FAC
        DYDX=DFAC
        XSQ=-0.5D0*X*X
        DXSQ=-X
        ISVAR=-(2*L+1)
        DO I=1,1000
          ISVAR=ISVAR+2
          DFAC=(DFAC*XSQ-X*FAC)/REAL(I*ISVAR,KIND=8)
          FAC=FAC*XSQ/REAL(I*ISVAR,KIND=8)
          Y=Y+FAC
          DYDX=DYDX+DFAC
          IF(ABS(FAC).LT.TOL) RETURN
        ENDDO
        CALL ERROR$MSG('Y NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE 
!       ========================================================================
!       ==  EXPANSION FOR LARGE ARGUMENTS   X>L                               ==
!       ========================================================================
        PI=4.D0*ATAN(1.D0)
        ARG=X+0.5D0*REAL(L,KIND=8)*PI
        M1POWERL=(-1)**L
        TRIG(1)=-M1POWERL*COS(ARG)/X
        DTRIG(1)=M1POWERL*(SIN(ARG)+COS(ARG)/X)/X
        TRIG(2)=M1POWERL*SIN(ARG)/X
        DTRIG(2)=M1POWERL*(COS(ARG)-SIN(ARG)/X)/X
        TRIG(3)=-TRIG(1)
        DTRIG(3)=-DTRIG(1)
        TRIG(4)=-TRIG(2)
        DTRIG(4)=-DTRIG(2)
        Y=TRIG(1)
        DYDX=DTRIG(1)
        IF(L.EQ.0) RETURN
!       ==  DOUBLE FACULTY FACUL(L)=(2*L)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*REAL(I,KIND=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
          DFAC=-FAC*REAL(K,KIND=8)/X
          Y=Y+FAC*TRIG(II)
          DYDX=DYDX+FAC*DTRIG(II)+DFAC*TRIG(II)
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$NEUMANN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODBESSEL(L,X,Y,dydx)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL BESSEL FUNCTION                   **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.2.5   FOR   X < L+1                                    **
!     **    FORMULA 10.2.9   FOR   X > L+1                                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dydx ! derivative
      REAL(8)               :: TRIG(2),dtrig(2)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ  !x-square
      REAL(8)               :: FAC,dfac
      INTEGER(4)            :: I,K,IL,ISVAR
!     **************************************************************************
      IF(X.GT.real(L+1,kind=8)) THEN
        TRIG(:)=1.D0
        dTRIG(:)=0.D0
        XSQ=0.5D0/X
        Y=XSQ*(EXP(X)-(-1.d0)**l*EXP(-X))
        dYdx=-y/x+xsq*(EXP(X)+(-1.d0)**l*EXP(-X))
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K,KIND=8)*XSQ
           dfac=-fac/x
           dtrig(1)=-dtrig(1)*fac-trig(1)*dfac
           TRIG(1)=-TRIG(1)*FAC
           dtrig(2)=dtrig(2)*fac+trig(2)*dfac
           TRIG(2)=+TRIG(2)*FAC
           dydx=dydx-XSQ/x*(TRIG(1)*EXP(X)-(-1)**l*TRIG(2)*EXP(-X)) &
    &                +XSQ*((dtrig(1)+TRIG(1))*EXP(X) &
               -(-1.d0)**l*(dtrig(2)-TRIG(2))*EXP(-X))
           Y=Y+XSQ*(TRIG(1)*EXP(X)-(-1.d0)**l*TRIG(2)*EXP(-X))
        ENDDO
!
!     ==========================================================================
!     ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                                ==
!     ==========================================================================
      ELSE
        ISVAR=1
        DO IL=1,L
          ISVAR=ISVAR*(2*IL+1)
        ENDDO
        if(l.eq.0) then
          FAC=1.D0/real(ISVAR,kind=8)
          dfac=0.d0
        else if(l.eq.1) then
          FAC=X/real(ISVAR,kind=8)
          dFAC=1.d0/real(ISVAR,kind=8)
        else
          FAC=X**L/real(ISVAR,kind=8)
          dFAC=real(l,kind=8)*X**(L-1)/real(ISVAR,kind=8)
        END IF
        Y=FAC
        dydx=dfac
        XSQ=0.5D0*X*X
        ISVAR=2*L+1
        DO I=1,1000
          ISVAR=ISVAR+2
          dFAC=(dFAC*XSQ+FAC*X)/real(I*ISVAR,kind=8)
          FAC=FAC*XSQ/real(I*ISVAR,kind=8)
          Y=Y+FAC
          dydx=dydx+dFAC
          IF(DABS(FAC).LT.TOL) RETURN
        ENDDO
        CALL ERROR$MSG('Y NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('spfunction$modbessel')
      END IF
      return
      END SUBROUTINE SPFUNCTION$MODBESSEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Spfunction$MODNEUMANN(L,X,Y,dydx)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL NEUMANN FUNCTION                  **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.2.6   FOR   X < L                                      **
!     **    FORMULA 10.2.10  FOR   X > L                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dydx ! derivative
      REAL(8)               :: TRIG(2),dtrig(2)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,dfac
      INTEGER(4)            :: I,K,IL,ISVAR
      logical(4)            :: convg
      real(8)               :: svar,dsvar
!     **************************************************************************
      if(x.le.0.d0) then
        CALL ERROR$MSG('not defined for zero or negative arguments')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('spfunction$modneumann')
      else if(X.le.real(L,kind=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                              ==
!       ========================================================================
        ISVAR= 1
        DO IL=2,L
          ISVAR=ISVAR*(2*IL-1)
        ENDDO
        FAC=real(ISVAR,kind=8)/X**(L+1)/(-1.d0)**L
        dfac=-real(l+1,kind=8)*fac/x
        Y=FAC
        dydx=dfac
        XSQ=0.5D0*X*X
        ISVAR=-(2*L+1)
        convg=.false.
        DO I=1,1000
          ISVAR=ISVAR+2
          dFAC=(dFAC*XSQ+FAC*X)/real(I*ISVAR,kind=8)
          FAC=FAC*XSQ/real(I*ISVAR,kind=8)
          Y=Y+FAC
          dydx=dydx+dFAC
          convg=ABS(FAC).LT.TOL
          if(convg) exit
        ENDDO
        if(.not.convg) then
          CALL ERROR$MSG('Y NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('X',X)
          CALL ERROR$STOP('spfunction$modneumann')
        end if
      else 
!       ========================================================================
!       ==  LARGE ARGUMENTS                                                   ==
!       ========================================================================
        TRIG(:)=1.D0
        DTRIG(:)=0.D0
        XSQ=0.5D0/X
        Y=XSQ*(EXP(X)+(-1.d0)**l*EXP(-X))
        DYDX=-Y/X+XSQ*(EXP(X)-(-1.d0)**l*EXP(-X))
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K, KIND=8)*XSQ
           DFAC=-FAC/X
           DTRIG(1)=-DTRIG(1)*FAC-TRIG(1)*DFAC
           TRIG(1)=-TRIG(1)*FAC
           DTRIG(2)= DTRIG(2)*FAC+TRIG(2)*DFAC
           TRIG(2)= TRIG(2)*FAC
           SVAR=XSQ*(TRIG(1)*DEXP(X)+(-1.D0)**l*TRIG(2)*DEXP(-X))
           DSVAR=-SVAR/X+XSQ*(           (DTRIG(1)+TRIG(1))*EXP(X) &
       &                     +(-1.D0)**l*(DTRIG(2)-TRIG(2))*EXP(-X))
           Y=Y+SVAR
           DYDX=DYDX+DSVAR
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODNEUMANN
!!$!
!!$!     .......................................................................
!!$      SUBROUTINE Spfunction$HANKEL(L,X,Y,dydx)
!!$!     ***********************************************************************
!!$!     **                                                                   **
!!$!     **  CALCULATES THE SPHERICAL HANKEL FUNCTION.                        **
!!$!     **  CONVENTION: ABRAMOWITZ, EQUATION 10.1.1                          **
!!$!     **                                                                   **
!!$!     ***********************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
!!$      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
!!$      COMPLEX(8),INTENT(OUT):: Y ! HANKEL FUNCTION AT X
!!$      COMPLEX(8),INTENT(OUT):: dydx ! HANKEL FUNCTION AT X
!!$      REAL(8)               :: YREAL, YIMAG,dYREAL, dYIMAG
!!$!     ***********************************************************************
!!$      CALL Spfunction$BESSEL (L,X,YREAL,dyreal)
!!$      CALL Spfunction$NEUMANN(L,X,YIMAG,dyimag)
!!$      Y    = CMPLX(YREAL,YIMAG)
!!$      dydx = CMPLX(dYREAL,dYIMAG)
!!$      END SUBROUTINE SPFUNCTION$HANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Spfunction$MODHANKEL(L,X,Y,dydx)
!     **************************************************************************
!     **  MODIFIED SPHERICAL BESSEL FUNCTION OF THE THIRD KIND K_L(X)         **
!     **  AS DEFINED IN ABRAMOWITZ/STEGUN EQ. 10.2.4                          **
!     **                                                                      **
!     **    Y = 0.5*PI*X^L*(-1/X D/DX)^L [EXP(-X)/X]                          **
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL HANKEL FUNCTION                   **
!     **  AS DEFINED BY ABRAMOWITZ AND STEGUN (EQUATION 10.2.4)               **
!     **  USING EQUATION 10.2.15 FOR X>L                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dYdx ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)               :: BESSEL_ABRA, NEUMANN_ABRA ! ABRAMOWITZ
      REAL(8)               :: PI
      REAL(8)               :: TRIG,dtrig
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,dfac
      REAL(8)               :: jval,jder
      REAL(8)               :: nval,nder
      REAL(8)               :: svar,dsvar
      INTEGER(4)            :: K
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      if(x.le.0.d0) then
        call error$msg('undefined for negative or zero arguments')
        call error$stop('Spfunction$MODHANKEL')
      else IF(X.le.real(L,kind=8)) THEN
         CALL SPFUNCTION$MODBESSEL (L,X,jval,jder)
         CALL SPFUNCTION$MODNEUMANN(L,X,nval,nder)
         Y = 0.5d0*PI * (-1.d0)**(L+1) * (jval - Nval)
         dydx = 0.5d0*PI * (-1.d0)**(L+1) * (jder - Nder)
      else 
        TRIG=1.D0
        dTRIG=0.D0
        XSQ=0.5D0/X
        Y=XSQ*PI*EXP(-X)
        dYdx=-y/x-y
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K,KIND=8)*XSQ
           dFAC=-fac/x
           dTRIG=dTRIG*FAC+TRIG*dFAC
           TRIG=TRIG*FAC
           svar=XSQ*PI*TRIG*EXP(-X)
           dsvar=-svar/x-svar+svar*dtrig/trig
           Y=Y+svar
           dYdx=dYdx+dsvar
        ENDDO
      END IF
      return
      END SUBROUTINE SPFUNCTION$MODHANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Spfunction$BESSEL0(L,X,Y,dydx)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION FOR K=0.                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dydx ! derivative
      INTEGER(4)            :: FAC, I
!     **************************************************************************
      IF(L.EQ.0) THEN
        Y = 1.D0
        DYDX=0.D0
      ELSE
        FAC=1.D0
        DO I=1, 2*L+1, 2
          FAC = FAC*I
        END DO
        Y = X**L / REAL(FAC,KIND=8)
        DYDX=REAL(L,KIND=8)*X**(L-1)/REAL(FAC,KIND=8)
      END IF
      END SUBROUTINE SPFUNCTION$BESSEL0
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Spfunction$NEUMANN0(L,X,Y,dydx)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION FOR K=0.                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: dYdx ! derivative
      INTEGER(4)            :: FAC, I
!     **************************************************************************
      if(x.le.0.d0) then
        call error$msg('not defined for zero or negative arguments')
        call error$stop('Spfunction$NEUMANN0')
      end if
      FAC=1.D0
      DO I=1,2*L-1,2
         FAC=FAC*I
      END DO
      Y = -real(FAC,kind=8) / X**(L+1)
      dydx = -real(l+1,kind=8)*y/x
      return
      END SUBROUTINE SPFUNCTION$NEUMANN0
