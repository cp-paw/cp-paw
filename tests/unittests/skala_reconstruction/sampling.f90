program skala_sampling
  use paw_skala_grid_module
  implicit none
  integer :: mode
  integer,parameter :: dimensions(3)=[1,2,4]
  do mode=1,size(dimensions)
    call check_sampling(dimensions(mode))
  end do
contains
subroutine check_sampling(nd)
  integer,intent(in) :: nd
  integer,parameter :: n1=9,n2=11,n3=13,np=7
  type(paw_skala_sample_cache) :: cache
  real(8) :: rho(n1,n2,n3,nd),grho(n1,n2,n3,3,nd),tau(n1,n2,n3,nd)
  real(8) :: cell(3,3),x(np,3),value(5,nd,np),jac(5,3,nd,np)
  real(8) :: r(nd),g(3,nd),t(nd),dr(3,nd),dg(3,3,nd),dt(3,nd),error
  real(8) :: backrho(n1,n2,n3,nd),backgrad(n1,n2,n3,3,nd),backtau(n1,n2,n3,nd)
  real(8) :: ar(np,max(2,nd)),ag(np,3,max(2,nd)),at(np,max(2,nd)),bar(5,nd)
  real(8) :: left,right,factor,h,xp(3),rp(nd),gp(3,nd),tp(nd),rm(nd),gm(3,nd),tm(nd)
  integer :: empty(0),index(1)
  integer :: i,j,k,d,c,p,trial,ng
  logical :: device,done,required,back_device
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
  ! One exactly coincident grid point also exercises the direct scatter branch.
  x(2,:)=matmul(cell,[2.d0/n1,3.d0/n2,4.d0/n3])
  index(1)=3+n1*(3+n2*4)
  do d=1,max(2,nd)
    do p=1,np
      ar(p,d)=sin(real(p+2*d,8))
      at(p,d)=cos(real(2*p-d,8))
      do c=1,3
        ag(p,c,d)=sin(real(p+2*c+3*d,8))
      end do
    end do
  end do
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
    left=0.d0
    do p=1,np
      do d=1,nd
        if(nd==4) then
          bar(1,d)=ar(p,d)
          bar(2:4,d)=ag(p,:,d)
          bar(5,d)=at(p,d)
        else
          factor=merge(1.d0,-1.d0,d==1)
          bar(1,d)=0.5d0*(ar(p,1)+factor*ar(p,2))
          bar(2:4,d)=0.5d0*(ag(p,:,1)+factor*ag(p,:,2))
          bar(5,d)=0.5d0*(at(p,1)+factor*at(p,2))
        end if
      end do
      left=left+sum(value(:,:,p)*bar)
    end do
    backrho=0.d0
    backgrad=0.d0
    backtau=0.d0
    call PAW_SKALA_GRID_BACK_BEGIN(n1,n2,n3,nd,backrho,backgrad,backtau,back_device)
    if(required.and..not.back_device) error stop 'required backprojection GPU was not used'
    call PAW_SKALA_GRID_BACKPROJECT(n1,n2,n3,cell,nd,0,np,empty,x,ar,ag,at, &
      backrho,backgrad,backtau,back_device)
    call PAW_SKALA_GRID_BACK_END(n1,n2,n3,nd,backrho,backgrad,backtau,back_device)
    right=sum(backrho*rho)+sum(backgrad*grho)+sum(backtau*tau)
    error=max(error,abs(left-right))
    if(trial==1) then
      ! Direct grid rows and interpolated rows must have the same transpose.
      backrho=0.d0
      backgrad=0.d0
      backtau=0.d0
      call PAW_SKALA_GRID_BACK_BEGIN(n1,n2,n3,nd,backrho,backgrad,backtau,back_device)
      call PAW_SKALA_GRID_BACKPROJECT(n1,n2,n3,cell,nd,1,1,index,x(2:2,:),ar(2:2,:), &
        ag(2:2,:,:),at(2:2,:),backrho,backgrad,backtau,back_device)
      call PAW_SKALA_GRID_BACK_END(n1,n2,n3,nd,backrho,backgrad,backtau,back_device)
      do d=1,nd
        factor=merge(1.d0,-1.d0,d==1)
        r(d)=ar(2,d)
        t(d)=at(2,d)
        g(:,d)=ag(2,:,d)
        if(nd/=4) then
          r(d)=0.5d0*(ar(2,1)+factor*ar(2,2))
          t(d)=0.5d0*(at(2,1)+factor*at(2,2))
          g(:,d)=0.5d0*(ag(2,:,1)+factor*ag(2,:,2))
        end if
      end do
      error=max(error,maxval(abs(backrho(3,4,5,:)-r)),maxval(abs(backtau(3,4,5,:)-t)), &
        maxval(abs(backgrad(3,4,5,:,:)-g)))
    end if
    h=1.d-5
    call PAW_SKALA_GRID_INTERPOLATE_PRIMITIVES(n1,n2,n3,cell,nd,rho,grho,tau,x(3,:),r,g,t,dr,dg,dt)
    do c=1,3
      xp=x(3,:)
      xp(c)=xp(c)+h
      call PAW_SKALA_GRID_INTERPOLATE_PRIMITIVES(n1,n2,n3,cell,nd,rho,grho,tau,xp,rp,gp,tp)
      xp(c)=xp(c)-2.d0*h
      call PAW_SKALA_GRID_INTERPOLATE_PRIMITIVES(n1,n2,n3,cell,nd,rho,grho,tau,xp,rm,gm,tm)
      if(maxval(abs((rp-rm)/(2.d0*h)-dr(c,:)))>2.d-8) error stop 'density coordinate derivative changed'
      if(maxval(abs((gp-gm)/(2.d0*h)-dg(:,c,:)))>2.d-8) error stop 'gradient coordinate derivative changed'
      if(maxval(abs((tp-tm)/(2.d0*h)-dt(c,:)))>2.d-8) error stop 'tau coordinate derivative changed'
    end do
    rho=0.9d0*rho
    grho=0.8d0*grho
    tau=0.7d0*tau
  end do
  call PAW_SKALA_SAMPLE_PREPARE(cache,n1,n2,n3,cell,x,0_8)
  if(PAW_SKALA_SAMPLE_BYTES(cache)/=0) error stop 'sample cache exceeded zero budget'
  call PAW_SKALA_GRID_SAMPLE_BATCH(cache,n1,n2,n3,nd,np,np,rho,grho,tau,value,jac,.false.,done)
  if(done) error stop 'disabled sample cache did not fall back'
  print *, 'SAMPLING VALUE/JACOBIAN/ADJOINT ERROR',nd,error
  if(error>1.d-10) error stop 'batched primitive interpolation changed'
  print *, 'PASS: sampling reuse, geometry invalidation, field updates and bounded fallback'
end subroutine check_sampling
end program skala_sampling
