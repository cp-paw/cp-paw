!
!     ..................................................................
      SUBROUTINE INBOUNDARY(THOM,AEZ,L,DEX,NR,R,U,UP,UPP,CF,CG,FR,FRP)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),   PARAMETER    :: FSS=(1.0D0/137.036D0)**2
      REAL(8),   INTENT(IN)   :: DEX
      REAL(8),   INTENT(IN)   :: AEZ
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: L ! MAIN ANGULAR MOMENTUM
      REAL(8),   INTENT(IN)   :: R(NR)
      REAL(8),   INTENT(INOUT):: U(NR)
      REAL(8),   INTENT(INOUT):: UP(NR)
      REAL(8),   INTENT(INOUT):: UPP(NR)
      REAL(8),   INTENT(IN)   :: CF(NR)
      REAL(8),   INTENT(IN)   :: CG(NR)
      REAL(8),   INTENT(IN)   :: FR(NR)
      REAL(8),   INTENT(IN)   :: FRP(NR)
      LOGICAL,  INTENT(IN)   :: THOM
      REAL(8)                 :: GAMMA
      REAL(8)                 :: RIGAMMA
      INTEGER(4)              :: IR, it !cwk
      REAL(8)                 :: U0(4),UP0(4), vtest !cwk
!     ******************************************************************
      IF(.NOT.THOM) THEN
         DO IR=1,4
            U(IR)=0.D0
            UP(IR)=0.D0
            UPP(IR)=CG(IR)
         ENDDO
         RETURN
      END IF
      IF(AEZ.GT.0) THEN
         IF(L .EQ. 0) THEN
            GAMMA=DSQRT(1.0D0-FSS*AEZ**2)
         ELSE
!$$$            GAMMA=(    L*DSQRT(   L**2 -FSS*AEZ**2)+ &
!$$$     &           (L+1)*DSQRT((L+1)**2-FSS*AEZ**2)   )/(2*L+1)
            gamma = sqrt(1.0+l*(l+1)-FSS*AEZ**2)
!cwk I have changed for gamma as defined in Koelling-Harmon
         END IF
      ELSE
         GAMMA=L+1
      END IF
! cwk your old machine
      
!!$      DO IR=1,4
!!$         RIGAMMA=R(ir)**GAMMA
!!$         U(IR)=RIGAMMA
!!$         UP(IR)=DEX*GAMMA*RIGAMMA
!!$         UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
!!$      ENDDO
!cwk my new  Attention!! U() means now u/r**gamma
      DO IR=1,4
        U(IR)=1.0
        UP(IR)=2.0*((gamma-1.0)/(fss*aez))*(gamma+2.0)/(2.0*gamma+1.0)*dex*r(ir)
!       print*,up(ir),-Aez*dex*r(ir) ! to jeszcze sprawdzic
        UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
     &          +(CF(IR)+FR(IR)+dex*gamma*(DEX+FRP(IR))-(dex*gamma)**2)*U(IR)+CG(IR)
      ENDDO

      u0=u(1:4)
      up0 = up(1:4)
      vtest = 1.0d2
      it=0
      do while (vtest .gt. 1.e-8)
        DO IR=2,4
          U(IR)=U(ir-1)+0.5*(UP(ir)+up(ir-1))
          UP(IR)=Up(ir-1)+0.5*(UPp(ir)+upp(ir-1))
          UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
    &            +(CF(IR)+FR(IR)+dex*gamma &
    &              *(DEX+FRP(IR))-(dex*gamma)**2)*U(IR)+CG(IR)
        ENDDO
        vtest=abs((u(4)-u0(4))/u0(4))+abs((up(4)-up0(4))/up0(4))
        u0=u(1:4)
        up0=up(1:4)
        it = it+1
        if (it .gt. 500 ) exit
      end do   
!     print *,'test in inbound', it, vtest !cwk
      do ir = 1,4
         RIGAMMA=R(ir)**GAMMA
!         print*,'bef',u(ir),up(ir),rigamma
!cwk here we are again at old U()
         up(ir)=(up(ir)+dex*u(ir)*gamma)*rigamma
         u(ir)=u(ir)*rigamma
   UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
!         print '(2x, 4f20.8)',u(ir),up(ir),upp(ir)
!!$         print *,(DEX+FRP(IR))*UP(IR),(CF(IR)+FR(IR))*U(IR)
      end do   
      RETURN

      END
