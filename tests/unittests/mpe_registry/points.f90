program skala_point_exchange
  use mpe_module
  use paw_skala_module
  implicit none
  integer, parameter :: np=7
  real(8) :: rho(np,2), grad(np,3,2), tau(np,2), x(np,3), w(np), d(np)
  real(8) :: weights(4)=[1.d0,4.d0,0.25d0,16.d0]
  integer :: nt, me, first, last, prior, lo, hi, n, k, a, b, c, mode
  integer, allocatable :: color(:)
  call MPE$INIT
  call MPE$QUERY('~',nt,me)
  allocate(color(nt))
  color=1
  call MPE$NEW('~','MONOMER',nt,color)
  point_tasks=nt
  point_rank=me
  point_cid='MONOMER'
  allocate(forward_deriv(3,1))
  do mode=1,2
    point_packed=mode==1
    do k=1,size(weights)
      point_root_weight=weights(k)
      do n=0,np
        prior=0
        do a=1,nt
          call SKALA_POINT_RANGE(n,a,lo,hi)
          if(lo/=prior+1.or.hi<lo-1) error stop 'noncontiguous point ranges'
          prior=hi
        end do
        if(prior/=n) error stop 'point range does not cover grid'
      end do
      call SKALA_POINT_RANGE(np,me,first,last)
      rho=0.d0
      grad=0.d0
      tau=0.d0
      x=0.d0
      w=0.d0
      d=0.d0
      do a=first,last
        do b=1,2
          rho(a,b)=100*b+a
          tau(a,b)=200*b+a
          do c=1,3
            grad(a,c,b)=1000*b+100*c+a
          end do
        end do
      end do
      call SKALA_POINT_EXCHANGE(np,.true.,rho,grad,tau,x,w,d)
      if(me==1) then
        call check(1,np)
        do a=1,np
          x(a,:)=a
          w(a)=2*a
          d(a)=3*a
        end do
      end if
      if(me/=1) then
        rho=0.d0
        grad=0.d0
        tau=0.d0
      end if
      call SKALA_POINT_EXCHANGE(np,.false.,rho,grad,tau,x,w,d)
      call check(first,last)
      do a=first,last
        if(any(x(a,:)/=a).or.w(a)/=2*a.or.d(a)/=3*a) error stop 'geometry scatter changed'
      end do
    end do
  end do
  if(me==1) print *, 'PASS: packed/unpacked exchange, weighted and empty point ranges'
  call MPE$EXIT
contains
  subroutine check(lo,hi)
    integer,intent(in) :: lo,hi
    integer :: i,j,l
    do i=lo,hi
      do j=1,2
        if(rho(i,j)/=100*j+i.or.tau(i,j)/=200*j+i) error stop 'field exchange changed'
        do l=1,3
          if(grad(i,l,j)/=1000*j+100*l+i) error stop 'gradient exchange changed'
        end do
      end do
    end do
  end subroutine check
end program skala_point_exchange
