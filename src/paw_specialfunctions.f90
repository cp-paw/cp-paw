!
!      .............................................................
       subroutine specialfunction__erf(x,val)
!      **  returns the error function erf(x)
!      ** see numerical recipes
       implicit none
       real(8),intent(in) :: x
       real(8),intent(out):: val
!      *************************************************************
       call specialfunction__gammp(x**2,0.5d0,val)
       If(x.lt.0.d0) val=-val
       return
       end
!
!      .............................................................
       subroutine specialfunction__gammp(x,a,val)
!      **  returns the incomplete gamma function P(a,x)
!      ** see numerical recipes
       implicit none
       real(8),intent(in) :: x
       real(8),intent(in) :: a
       real(8),intent(out):: val
!      *************************************************************
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
!      .............................................................
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
!      .............................................................
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
!      .............................................................
       SUBROUTINE SPECIALFUNCTION_Gcf(X,A,val)
!      ** 
!      ** incomplete gamma function q(a,x) evaluated by its continued fraction representation
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
!     ..............................................BESSL...............
      SUBROUTINE SPECIALFUNCTION$BESSEL(L,X,Y)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                    **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                          **
!     **    FORMULA 10.1.2              FOR   X < 8                   **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > 8                   **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
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
      IF(X.GT.DBLE(L)) THEN
        PI=4.D0*DATAN(1.D0)
        ARG=X-0.5D0*DBLE(L)*PI
        TRIG(1)=DSIN(ARG)/X
        TRIG(2)=DCOS(ARG)/X
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        Y=TRIG(1)
        IF(L.EQ.0) RETURN
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*DBLE(I)
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
        IF(DABS(FAC).LT.TOL) GOTO 9999
      ENDDO
      CALL ERROR$MSG('Y NOT CONVERGED')
      CALL ERROR$I4VAL('L',L)
      CALL ERROR$R8VAL('X',X)
      CALL ERROR$STOP('specialfunction$bessel')
9999  CONTINUE
      RETURN
      END
