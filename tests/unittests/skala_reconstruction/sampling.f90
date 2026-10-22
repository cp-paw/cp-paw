program skala_sampling
  use paw_skala_grid_module
  implicit none
  integer,parameter :: n1=9,n2=11,n3=13,nd=2,np=7
  type(paw_skala_sample_cache) :: cache
  real(8) :: rho(n1,n2,n3,nd),grho(n1,n2,n3,3,nd),tau(n1,n2,n3,nd)
  real(8) :: cell(3,3),x(np,3),value(5,nd,np),jac(5,3,nd,np)
  real(8) :: r(nd),g(3,nd),t(nd),dr(3,nd),dg(3,3,nd),dt(3,nd),error
  integer :: i,j,k,d,c,p,trial,ng
  logical :: device,done,required
  required=command_argument_count()>0
  cell=reshape([4.d0,0.d0,0.d0,0.3d0,5.d0,0.d0,-0.2d0,0.1d0,6.d0],[3,3])
  do d=1,nd
    do k=1,n3
      do j=1,n2
        do i=1,n1
          rho(i,j,k,d)=sin(real(i+2*j+3*k,8))*d
          tau(i,j,k,d)=cos(real(2*i-j+k,8))*d
          do c=1,3
            grho(i,j,k,c,d)=sin(real(c*i-j+2*k,8))*d
          end do
        end do
      end do
    end do
  end do
  do p=1,np
    x(p,:)=real(p,8)*[0.2d0,-0.7d0,1.1d0]
  end do
  x(1,:)=[-1.d-12,5.d0,12.d0]
  error=0.d0
  do trial=1,4
    if(trial==3) x(2,1)=x(2,1)+0.01d0
    if(trial==4) cell(1,2)=cell(1,2)+0.02d0
    call PAW_SKALA_SAMPLE_PREPARE(cache,n1,n2,n3,cell,x,1000000_8)
    if(PAW_SKALA_SAMPLE_BYTES(cache)<=0) error stop 'sample cache did not allocate'
    ng=np
    if(trial==2) ng=0
    call PAW_SKALA_GRID_SAMPLE_BEGIN(n1,n2,n3,nd,rho,grho,tau,device)
    if(required.and..not.device) error stop 'required sampling GPU was not used'
    call PAW_SKALA_GRID_SAMPLE_BATCH(cache,n1,n2,n3,nd,np,ng,rho,grho,tau,value,jac(:,:,:,1:ng),device,done)
    if(.not.done) error stop 'sample batch fell back unexpectedly'
    call PAW_SKALA_GRID_SAMPLE_END(n1,n2,n3,nd,rho,grho,tau,device)
    do p=1,np
      call PAW_SKALA_GRID_INTERPOLATE_PRIMITIVES(n1,n2,n3,cell,nd,rho,grho,tau,x(p,:),r,g,t,dr,dg,dt)
      error=max(error,maxval(abs(value(1,:,p)-r)),maxval(abs(value(2:4,:,p)-g)),maxval(abs(value(5,:,p)-t)))
      if(ng>0) error=max(error,maxval(abs(jac(1,:,:,p)-dr)), &
        maxval(abs(jac(2:4,:,:,p)-dg)),maxval(abs(jac(5,:,:,p)-dt)))
    end do
    rho=0.9d0*rho
    grho=0.8d0*grho
    tau=0.7d0*tau
  end do
  call PAW_SKALA_SAMPLE_PREPARE(cache,n1,n2,n3,cell,x,0_8)
  if(PAW_SKALA_SAMPLE_BYTES(cache)/=0) error stop 'sample cache exceeded zero budget'
  call PAW_SKALA_GRID_SAMPLE_BATCH(cache,n1,n2,n3,nd,np,np,rho,grho,tau,value,jac,.false.,done)
  if(done) error stop 'disabled sample cache did not fall back'
  print *, 'SAMPLING VALUE/JACOBIAN ERROR',error
  if(error>1.d-10) error stop 'batched primitive interpolation changed'
  print *, 'PASS: sampling reuse, geometry invalidation, field updates and bounded fallback'
end program skala_sampling
