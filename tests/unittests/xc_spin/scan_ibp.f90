program scan_ibp
  use paw_spin_quadrature_module, only: paw_spin_polar_ibp_grid
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(8), allocatable :: axes(:,:),qw(:),coef(:),resp(:,:)
  real(8), allocatable :: folded_axes(:,:),folded_w(:),folded_c(:),folded_b(:,:)
  real(8) :: x(21),y(21),d(21),tmp(21),axis(3),v(3),q,e,ep,em,h,err,fd
  integer :: i,k
  call MPE$INIT
  call SCAN_INITIALIZE()
  call paw_spin_polar_ibp_grid(16,16,axes,qw,coef,resp)
  call paw_spin_polar_ibp_grid(16,16,folded_axes,folded_w,folded_c,folded_b,antipodal_even=.true.)
  if(size(folded_w)/=257) error stop 'Folded IBP grid count'
  if(size(qw)/=514) error stop 'IBP grid count'
  v=[.3d0,-.7d0,.2d0]
  do k=1,4
    e=0.d0
    do i=1,size(qw)
      q=dot_product(v,axes(:,i))
      e=e+qw(i)*(coef(i)*q**k+k*q**(k-1)*dot_product(v,resp(:,i)))
    end do
    if(abs(e-merge(sum(v**2)**(k/2),0.d0,mod(k,2)==0))>2.d-12) error stop 'IBP angular identity'
  end do
  axis=[1.d0,2.d0,-3.d0]/sqrt(14.d0)
  x=[.3d0,.24d0*axis,.03d0,.02d0,-.01d0, &
     .8d0*axis(1)*[.03d0,.02d0,-.01d0],.8d0*axis(2)*[.03d0,.02d0,-.01d0], &
     .8d0*axis(3)*[.03d0,.02d0,-.01d0],.2d0,.16d0*axis,1.d0]
  call eval(x,e,d)
  ! Independent Libxc SCAN exchange and correlation, FHC disabled.
  if(abs(e-(-.1979705104499099d0))>1.d-10) error stop 'Rotated collinear SCAN reference'
  x(8:10)=x(8:10)+[.017d0,-.01d0,.012d0]
  x(18:20)=x(18:20)+[.013d0,.021d0,.001d0]
  x(21)=.69d0
  call eval(x,e,d)
  err=0.d0
  h=1.d-6
  do i=1,21
    y=x
    y(i)=x(i)+h
    call eval(y,ep,tmp)
    y(i)=x(i)-h
    call eval(y,em,tmp)
    fd=(ep-em)/(2.d0*h)
    if(.not.ieee_is_finite(fd).or..not.ieee_is_finite(d(i))) error stop 'Nonfinite SCAN IBP response'
    err=max(err,abs(fd-d(i)))
  end do
  print '(A,ES14.5)','SCAN POLAR_IBP MAX FINITE-DIFFERENCE ERROR ',err
  if(err>2.d-8) error stop 'SCAN IBP field or weight derivative'
  call MPE$EXIT
contains
  subroutine eval(z,en,der)
    real(8),intent(in) :: z(21)
    real(8),intent(out) :: en,der(21)
    real(8) :: dg(3,4),fg(3,4),fe,fr(4),ft(4),fw,err_fold
    call scan_polar_point(z(1:4),reshape(z(5:16),[3,4]),z(17:20),z(21),size(qw), &
      axes,qw,coef,resp,en,der(1:4),dg,der(17:20),der(21))
    der(5:16)=reshape(dg,[12])
    call scan_polar_point(z(1:4),reshape(z(5:16),[3,4]),z(17:20),z(21),size(folded_w), &
      folded_axes,folded_w,folded_c,folded_b,fe,fr,fg,ft,fw)
    err_fold=max(abs(en-fe),maxval(abs(der-[fr,reshape(fg,[12]),ft,fw])))
    if(.not.ieee_is_finite(err_fold).or.err_fold>2.d-12) error stop 'Antipodal energy or pullback'
  end subroutine
end program
