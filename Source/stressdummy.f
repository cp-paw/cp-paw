!
!     ..................................................................
      subroutine potential_stresstest
      implicit none
      integer(4),parameter :: lmrxx=9
      real(8)              :: a(3,3,lmrxx,lmrxx)
      real(8)              :: a1(3,3)
      real(8)              :: fourpiby3
      real(8)              :: dsmall=1.d+12
      real(8)              :: cg1,cg2
      integer(4),parameter :: maplm(3)=(/2,4,3/)
!     ******************************************************************
      fourpiby3=4.d0*(4.d0*datan(1.d0)/3.d0
      do lm1=1,lmrxx
        l1=int(sqrt(lm1-1)+dsmall)
        lm3x=(l1+1)**2
        do lm2=1,lmrxx
          a1(:,:)=0.d0
          do lm3=1,lm3x
            do i=1,3
              lmi=maplm(i)
              do i=1,3
                lmj=maplm(i)
                call clebsch(lmi,lm1,lm3,cg1)
                call clebsch(lmj,lm2,lm3,cg2)
                a1(i,j)=a1(i,j)+cg1*cg2
              enddo
            enddo
          enddo 
          a(:,:,lm1,lm2)=fourpiby3*real(2*l1+1,kind=8)*a1(:,:)
        enddo
      enddo
      DO I=1,3
        DO J=1,3
          DO LM1=1,LMRXX
            WRITE(*,FMT='(9F10.5)')A(I,J,LM1,:)
          ENDDO
        ENDDO
      ENDDO
      return
      end
