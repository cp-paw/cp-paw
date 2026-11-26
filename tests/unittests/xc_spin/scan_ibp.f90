program scan_ibp
  use paw_spin_quadrature_module, only: paw_spin_polar_ibp_grid,paw_spin_frame_prepare,paw_spin_frame
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(8), allocatable :: axes(:,:),qw(:),coef(:),resp(:,:)
  real(8), allocatable :: folded_axes(:,:),folded_w(:),folded_c(:),folded_b(:,:)
  real(8) :: x(21),y(21),d(21),tmp(21),axis(3),v(3),q,e,ep,rot(3,3),want(21),aligned(21)
  type(paw_spin_frame) :: frame
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
  aligned=x
  x(8:10)=x(8:10)+[.017d0,-.01d0,.012d0]
  x(18:20)=x(18:20)+[.013d0,.021d0,.001d0]
  x(21)=.69d0
  call eval(x,e,d)
  ! Rotate only Pauli components. Spatial gradient axes and the weight stay fixed.
  rot=reshape([-20, -4, 22, -20, -10, -20, 10, -28, 4],[3,3])/30.d0
  call rotate(x,y)
  call rotate(d,want)
  call eval(y,ep,tmp)
  if(abs(e-ep)>2.d-12.or.maxval(abs(want-tmp))>2.d-10) error stop 'Noncollinear SCAN spin covariance'
  call check_response(x,'resolved')
  aligned(21)=x(21)
  x=aligned+.15d0*(x-aligned)
  call paw_spin_frame_prepare(1,reshape(x(1:4),[1,4]),reshape(x(5:16),[1,3,4]), &
                              reshape(x(17:20),[1,4]),[1.d0],frame)
  if(.not.(frame%gauge_weight(3)>.25d0.and.frame%gauge_weight(3)<.75d0)) error stop 'Expected frame transition'
  call check_response(x,'transition')
  call MPE$EXIT
contains
  subroutine check_response(z,label)
    real(8),intent(in) :: z(21)
    character(*),intent(in) :: label
    real(8),parameter :: h=1.d-6
    real(8) :: zz(21),response(21),other(21),en,plus,minus,err
    integer :: j
    call eval(z,en,response)
    err=0.d0
    do j=1,21
      zz=z
      zz(j)=z(j)+h
      call eval(zz,plus,other)
      zz(j)=z(j)-h
      call eval(zz,minus,other)
      err=max(err,abs((plus-minus)/(2.d0*h)-response(j)))
    end do
    print '(3A,ES14.5)','SCAN POLAR_IBP ',label,' MAX FINITE-DIFFERENCE ERROR ',err
    if(.not.ieee_is_finite(err).or.err>2.d-8) error stop 'SCAN IBP field or weight derivative'
  end subroutine
  subroutine rotate(a,b)
    real(8),intent(in) :: a(21)
    real(8),intent(out) :: b(21)
    integer :: j
    b=a
    b(2:4)=matmul(rot,a(2:4))
    b(18:20)=matmul(rot,a(18:20))
    do j=8,10
      b(j:j+6:3)=matmul(rot,a(j:j+6:3))
    end do
  end subroutine
  subroutine eval(z,en,der)
    real(8),intent(in) :: z(21)
    real(8),intent(out) :: en,der(21)
    real(8) :: dg(3,4),fg(3,4),fe,fr(4),ft(4),fw,err_fold
    call scan_polar_point(z(1:4),reshape(z(5:16),[3,4]),z(17:20),z(21),size(qw), &
      axes,qw,coef,resp,en,der(1:4),dg,der(17:20),der(21))
    der(5:16)=reshape(dg,[12])
    call scan_polar_point(z(1:4),reshape(z(5:16),[3,4]),z(17:20),z(21),size(folded_w), &
      folded_axes,folded_w,folded_c,folded_b,fe,fr,fg,ft,fw)
    if(.not.all(ieee_is_finite([en,der,fe,fr,reshape(fg,[12]),ft,fw]))) error stop 'Nonfinite SCAN IBP response'
    err_fold=max(abs(en-fe),maxval(abs(der-[fr,reshape(fg,[12]),ft,fw])))
    if(.not.ieee_is_finite(err_fold).or.err_fold>2.d-12) error stop 'Antipodal energy or pullback'
  end subroutine
end program
