program grid
  use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  interface
    integer(c_int) function setenv(name,value,overwrite) bind(C,name="setenv")
      import c_int,c_char
      character(c_char), intent(in) :: name(*),value(*)
      integer(c_int),value :: overwrite
    end function
  end interface
  integer, parameter :: n=4103,nx=n+3
  real(8) :: rho(nx,2),g(n,3,2),vh(nx,2),vg(n,3,2),v(nx,2),dg(n,3,2)
  real(8) :: eh,ed,sh(3,3),sd(3,3),delta(4),d1(5),d2(5,5),d3(5,5,5)
  integer :: i,s,k,rc
  character(8) :: family
  logical :: active,tgra
  call get_command_argument(1,family)
  tgra=family=='PBE'
  if (.not.tgra.and.family/='VWN') error stop 'select PBE or VWN'
  call DFT$SETI4('TYPE',merge(111111,20,tgra))
  rho=0.d0
  do i=1,n
    rho(i,1)=0.4d0+0.15d0*sin(0.03d0*i)
    rho(i,2)=rho(i,1)*0.45d0*cos(0.07d0*i)
    do k=1,3
      g(i,k,1)=0.04d0*cos(0.02d0*i+k)
      g(i,k,2)=0.02d0*sin(0.03d0*i+k)
    end do
  end do
  rho(1,1)=-1.d-9
  rho(2,1)=4.d-7
  rho(n,2)=2.d0*rho(n,1)
  g(1:2,:,:)=0.d0
  do s=1,2
    call mode('cpu')
    call PAWLIBXC$CUDA_ACTIVE(active)
    if(active) error stop 'host reference executed on device'
    vh=rho
    vg=g
    call POTENTIAL_XC(tgra,s,nx,n,n,17.d0,vh,vg,eh,.true.,sh)
    call mode('cuda')
    call PAWLIBXC$CUDA_ACTIVE(active)
    if(.not.active) error stop 'device test silently fell back to host'
    v=rho
    dg=g
    call POTENTIAL_XC(tgra,s,nx,n,n,17.d0,v,dg,ed,.true.,sd)
    delta=[abs(eh-ed),maxval(abs(v(:,1:s)-vh(:,1:s))), &
           maxval(abs(dg(:,:,1:s)-vg(:,:,1:s))),maxval(abs(sh-sd))]
    print '(a,i2,4es14.5)',trim(family),s,delta
    if(.not.all(ieee_is_finite(delta)).or.any(delta>2.d-10)) error stop 'host/device XC mismatch'
  end do
  call PAWLIBXC$GGA3([0.5d0,0.12d0,0.05d0,0.01d0,0.01d0],eh,d1,d2,d3)
  if(.not.all(ieee_is_finite(d3)).or.maxval(abs(d3))==0.d0) error stop 'host KXC unavailable'
  rc=setenv('CPPAW_GPU_MODE'//c_null_char,'off'//c_null_char,1_c_int)
  if(rc/=0) error stop 'setenv failed'
  call PAWLIBXC$CUDA_ACTIVE(active)
  if(active) error stop 'GPU off did not disable Libxc'
contains
  subroutine mode(value)
    character(*),intent(in) :: value
    rc=setenv('CPPAW_LIBXC_DEVICE'//c_null_char,value//c_null_char,1_c_int)
    if(rc/=0) error stop 'setenv failed'
  end subroutine
end program
