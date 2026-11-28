program xc_spin_noncollinear
  use paw_lebedev_module, only: paw_lebedev_size,paw_lebedev_get
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: np=16,na=2
  real(8) :: r(np,4),g(np,3,4),t(np,4),r0(np,4),g0(np,3,4),t0(np,4)
  real(8) :: dr(np,4),dg(np,3,4),dt(np,4),dr0(np,4),dg0(np,3,4),dt0(np,4)
  real(8) :: x(np,3),w(np),c(na,3),aw(np),dx(np,3),dw(np),dc(na,3),daw(np)
  real(8) :: x0(np,3),w0(np),c0(na,3),aw0(np),dx0(np,3),dw0(np),dc0(na,3),daw0(np)
  real(8) :: ar(np,4),ag(np,3,4),at(np,4),e,e0,ep,em,h,fd,analytic,error,tolerance
  real(8) :: rot(3,3),axis(3),angle,radius,steps(2),rc(np,2),gc(np,3,2),tc(np,2),ec
  real(8) :: drc(np,2),dgc(np,3,2),dtc(np,2),bound
  real(8) :: dxc(np,3),dwc(np),dcc(na,3),dawc(np)
  real(8) :: stats(8)
  real(8),allocatable :: axes(:,:),weights(:),original_axes(:,:)
  integer(8) :: sizes(na)
  integer :: i,j,k,s,nq,ia,level
  character(16) :: functional,precision,device
  character(2048) :: model
  logical :: available
  if(command_argument_count()<1.or.command_argument_count()>4) &
    error stop 'Usage: xc_noncollinear.x SCAN|SKALA [MODEL [published|diagnostic64 [CPU|CUDA]]]'
  call get_command_argument(1,functional)
  if(functional/='SCAN'.and.functional/='SKALA') error stop 'Unknown XC functional'
  model=''
  precision='published'
  device='CPU'
  if(command_argument_count()>=2) call get_command_argument(2,model)
  if(command_argument_count()>=3) call get_command_argument(3,precision)
  if(command_argument_count()>=4) call get_command_argument(4,device)
  if(functional=='SKALA'.and.len_trim(model)==0) error stop 'Skala requires a model'
  if(precision/='published'.and.precision/='diagnostic64') error stop 'Unknown model precision'
  call MPE$INIT
  call SKALA$SETCH('FUNCTIONAL',trim(functional))
  call SKALA$SETL4('ON',.true.)
  call SKALA$SETL4('APPLY',.true.)
  call SKALA$SETL4('NEEDFORCE',.true.)
  if(functional=='SKALA') then
    call SKALA$SETCH('MODEL',trim(model))
    call SKALA$SETCH('DEVICE',trim(device))
  end if
  call SKALA$AVAILABLE(available)
  if(.not.available) error stop 'Requested XC backend is unavailable'
  c=0.d0
  c(2,3)=1.4d0
  sizes=np/na
  do i=1,np
    ia=1+(i-1)/(np/na)
    j=1+mod(i-1,np/na)
    angle=2.d0*acos(-1.d0)*real(j-1,8)/real(np/na,8)
    radius=0.35d0+0.03d0*mod(j,3)
    x(i,:)=c(ia,:)+[radius*cos(angle),radius*sin(angle),0.02d0*(j-4)]
    r(i,1)=0.5d0*exp(-radius)
    r(i,2:4)=0.1d0*r(i,1)*[cos(angle),sin(angle),0.5d0]
    g(i,:,1)=[-0.24d0*cos(angle),-0.24d0*sin(angle),0.02d0*sin(angle)]
    do k=2,4
      g(i,:,k)=0.013d0*[sin(angle+k),cos(angle-k),sin(2.d0*angle+k)]
    end do
    t(i,1)=0.18d0*r(i,1)
    t(i,2:4)=0.07d0*t(i,1)*[sin(angle),0.4d0,cos(angle)]
    w(i)=0.04d0+0.001d0*j
    aw(i)=0.05d0
  end do
  r0=r
  g0=g
  t0=t
  x0=x
  w0=w
  c0=c
  aw0=aw
  call quadrature(3)
  call evaluate
  if(functional=='SKALA') then
    ! TorchScript may specialize its graph during the first evaluations.
    ! Retain and bound the first-call difference before the strict dispatch
    ! comparison between warmed calls with otherwise identical inputs.
    e0=e
    dr0=dr
    dg0=dg
    dt0=dt
    do level=1,3
      call evaluate
    end do
    error=max(abs(e-e0),maxval(abs(dr-dr0)),maxval(abs(dg-dg0)),maxval(abs(dt-dt0)))
    write(*,'(A,ES18.10)') 'PAULI XC FIRST-CALL TO WARMED DIFFERENCE ',error
    bound=5.d-8
    if(precision=='published') bound=2.d-6
    if(error>bound) error stop 'Skala graph warmup changed the result beyond its precision bound'
  end if
  if(functional=='SCAN') then
    call SCAN$PAULI_FIELD_LIMITS(np,nq,r,g,t,w,axes,weights,stats)
    if(abs(stats(1)-2.d0*np)>1.d-12.or.abs(stats(6)-sum(w*r(:,1)))>1.d-12) &
      error stop 'Projected SCAN statistics count spin directions as extra electrons'
  end if
  e0=e
  dr0=dr
  dg0=dg
  dt0=dt
  call SKALA$SETI4('SPINEXACTNESS',3)
  call SKALA_EVALUATE_FIELDS(np,na,4,r,g,t,x,w,c,aw,sizes,e,dr,dg,dt,dx,dw,dc,daw)
  write(*,'(A,4ES18.10)') 'PAULI XC DISPATCH ENERGY/RHO/GRAD/TAU ERRORS ', &
    abs(e-e0),maxval(abs(dr-dr0)),maxval(abs(dg-dg0)),maxval(abs(dt-dt0))
  if(max(abs(e-e0),maxval(abs(dr-dr0)),maxval(abs(dg-dg0)),maxval(abs(dt-dt0)))>1.d-12) &
    error stop 'PAW XC dispatch differs from explicit Pauli backend'
  e0=e
  dr0=dr
  dg0=dg
  dt0=dt
  dx0=dx
  dw0=dw
  dc0=dc
  daw0=daw
  if(abs(e)<1.d-12.or.maxval(abs(dr(:,2:4)))<1.d-12) error stop 'Degenerate Pauli XC result'
  steps=[1.d-4,3.d-5]
  tolerance=5.d-8
  if(functional=='SCAN') tolerance=5.d-7
  if(functional=='SKALA'.and.precision=='published') then
    steps=[1.d-3,3.d-4]
    tolerance=2.d-5
  end if
  do k=1,4
    ar=0.d0
    ag=0.d0
    at=0.d0
    if(k==1) then
      ar(:,1)=0.07d0*r0(:,1)
      ag(:,:,1)=-0.04d0*g0(:,:,1)
      at(:,1)=0.09d0*t0(:,1)
    else if(k==2) then
      do s=2,4
        ar(:,s)=0.03d0*r0(:,1)*sin(real(s,8))
        ag(:,:,s)=0.04d0*g0(:,:,1)*cos(real(s,8))
      end do
    else if(k==3) then
      do s=2,4
        at(:,s)=0.11d0*t0(:,1)*sin(real(s,8))
      end do
    end if
    analytic=sum(dr0*ar)+sum(dg0*ag)+sum(dt0*at)
    if(k==4) analytic=0.03d0*(sum(dx0*x0)+sum(dw0*w0)+sum(dc0*c0)+sum(daw0*aw0))
    do j=1,size(steps)
      h=steps(j)
      call perturb(h)
      call evaluate
      ep=e
      call perturb(-h)
      call evaluate
      em=e
      fd=(ep-em)/(2.d0*h)
      error=abs(fd-analytic)/max(1.d0,abs(analytic))
      write(*,'(A,I2,A,ES12.4,A,ES18.10)') 'PAULI XC DIRECTION ',k,' STEP ',h,' FD ERROR ',error
      if(error>tolerance) error stop 'Multicollinear XC adjoint mismatch'
    end do
  end do
  call perturb(0.d0)
  ! Rotate spin fields and quadrature together. This covariance is exact even
  ! for a finite spin quadrature and tests every returned Pauli component.
  original_axes=axes
  axis=[1.d0,2.d0,-3.d0]/sqrt(14.d0)
  angle=0.71d0
  rot=0.d0
  do i=1,3
    rot(i,i)=cos(angle)
    do j=1,3
      rot(i,j)=rot(i,j)+(1.d0-cos(angle))*axis(i)*axis(j)
    end do
  end do
  rot(1,2)=rot(1,2)-sin(angle)*axis(3)
  rot(2,1)=rot(2,1)+sin(angle)*axis(3)
  rot(1,3)=rot(1,3)+sin(angle)*axis(2)
  rot(3,1)=rot(3,1)-sin(angle)*axis(2)
  rot(2,3)=rot(2,3)-sin(angle)*axis(1)
  rot(3,2)=rot(3,2)+sin(angle)*axis(1)
  r(:,2:4)=matmul(r0(:,2:4),transpose(rot))
  t(:,2:4)=matmul(t0(:,2:4),transpose(rot))
  do k=1,3
    g(:,k,2:4)=matmul(g0(:,k,2:4),transpose(rot))
  end do
  axes=matmul(rot,original_axes)
  call evaluate
  error=max(abs(e-e0),maxval(abs(dr(:,1)-dr0(:,1))),maxval(abs(dt(:,1)-dt0(:,1))), &
    maxval(abs(dg(:,:,1)-dg0(:,:,1))),maxval(abs(dr(:,2:4)-matmul(dr0(:,2:4),transpose(rot)))), &
    maxval(abs(dt(:,2:4)-matmul(dt0(:,2:4),transpose(rot)))),maxval(abs(dx-dx0)), &
    maxval(abs(dw-dw0)),maxval(abs(dc-dc0)),maxval(abs(daw-daw0)))
  do k=1,3
    error=max(error,maxval(abs(dg(:,k,2:4)-matmul(dg0(:,k,2:4),transpose(rot)))))
  end do
  write(*,'(A,ES18.10)') 'PAULI XC ROTATION COVARIANCE ERROR ',error
  bound=5.d-8
  if(functional=='SKALA'.and.precision=='published') bound=2.d-6
  if(error>bound) error stop 'Multicollinear spin-rotation covariance'
  call perturb(0.d0)
  r(:,2:4)=0.d0
  call evaluate
  if(maxval(abs(dt(:,2:4)))<1.d-12) error stop 'Spin kinetic density was lost at zero magnetization'
  g(:,:,2:4)=0.d0
  t(:,2:4)=0.d0
  do s=1,2
    rc(:,s)=0.5d0*r(:,1)
    gc(:,:,s)=0.5d0*g(:,:,1)
    tc(:,s)=0.5d0*t(:,1)
  end do
  call SKALA$EVALUATE(np,na,rc,gc,tc,x,w,c,aw,sizes,ec,drc,dgc,dtc,dxc,dwc,dcc,dawc)
  call evaluate
  error=max(abs(e-ec),maxval(abs(dr(:,1)-drc(:,1))),maxval(abs(dg(:,:,1)-dgc(:,:,1))), &
    maxval(abs(dt(:,1)-dtc(:,1))),maxval(abs(dr(:,2:4))),maxval(abs(dg(:,:,2:4))), &
    maxval(abs(dt(:,2:4))),maxval(abs(dx-dxc)),maxval(abs(dw-dwc)), &
    maxval(abs(dc-dcc)),maxval(abs(daw-dawc)))
  write(*,'(A,ES18.10)') 'PAULI XC CLOSED-SHELL LIMIT ERROR ',error
  bound=5.d-8
  if(functional=='SKALA'.and.precision=='published') bound=2.d-6
  if(error>bound) error stop 'Multicollinear closed-shell limit'
  ! The collinear limit is an angular-convergence test, not the six-axis rule.
  r=r0
  g=g0
  t=t0
  r(:,2:3)=0.d0
  g(:,:,2:3)=0.d0
  t(:,2:3)=0.d0
  do s=1,2
    rc(:,s)=0.5d0*(r(:,1)+real(3-2*s,8)*r(:,4))
    gc(:,:,s)=0.5d0*(g(:,:,1)+real(3-2*s,8)*g(:,:,4))
    tc(:,s)=0.5d0*(t(:,1)+real(3-2*s,8)*t(:,4))
  end do
  call SKALA$EVALUATE(np,na,rc,gc,tc,x,w,c,aw,sizes,ec,drc,dgc,dtc,dxc,dwc,dcc,dawc)
  do level=17,53,36
    call quadrature(level)
    call evaluate
    error=max(abs(e-ec),maxval(abs(dr(:,1)-0.5d0*(drc(:,1)+drc(:,2)))), &
      maxval(abs(dr(:,4)-0.5d0*(drc(:,1)-drc(:,2)))), &
      maxval(abs(dg(:,:,1)-0.5d0*(dgc(:,:,1)+dgc(:,:,2)))), &
      maxval(abs(dg(:,:,4)-0.5d0*(dgc(:,:,1)-dgc(:,:,2)))), &
      maxval(abs(dt(:,1)-0.5d0*(dtc(:,1)+dtc(:,2)))), &
      maxval(abs(dt(:,4)-0.5d0*(dtc(:,1)-dtc(:,2)))),maxval(abs(dr(:,2:3))), &
      maxval(abs(dg(:,:,2:3))),maxval(abs(dt(:,2:3))),maxval(abs(dx-dxc)), &
      maxval(abs(dw-dwc)),maxval(abs(dc-dcc)),maxval(abs(daw-dawc)))
    write(*,'(A,I4,A,ES18.10)') 'PAULI XC COLLINEAR LIMIT ',nq,' ERROR ',error
    bound=5.d-7
    if(functional=='SKALA'.and.precision=='published') bound=2.d-6
    if(error>bound) error stop 'Multicollinear collinear limit'
  end do
  if(functional=='SCAN') then
    call vacuum_response
    call strong_scan_reference
  end if
  call SKALA$FINALIZE
  write(*,'(A)') trim(functional)//' MULTICOLLINEAR BACKEND CHECKS PASSED'
  call MPE$EXIT
contains

  subroutine strong_scan_reference
    ! Frozen independent MCfun 0.2.6 / LibXC 7.0.0 SCAN 263+267 values.
    ! LibXC Fermi-hole-curvature enforcement and MCfun's local-collinear
    ! fallback were disabled. Identical 110-node Lebedev quadratures are used.
    ! The pole-aligned p=.9999 case is ill-conditioned in double precision.
    real(8),parameter :: polarization(3)=[0.8d0,0.9999d0,0.9999d0]
    real(8),parameter :: eref(3)=[-1.98023028619234709d-1, -2.13405668891177791d-1, -2.13256966851859242d-1]
    real(8),parameter :: reference(5,4,3)=reshape([ &
      -7.38067623213678914d-1, -1.12005224467018262d-3, -7.46701496446788631d-4, &
      3.73350748223394424d-4, 3.20399735116302931d-2, &
      0.00000000000000000d+0, 0.00000000000000000d+0, 0.00000000000000000d+0, &
      0.00000000000000000d+0, 0.00000000000000000d+0, &
      4.33680868994201774d-19, -2.03287907341032081d-20, -6.77626357803440271d-21, &
      3.38813178901720136d-21, -6.50521303491302660d-19, &
      -2.17847969435823074d-1, 4.34718288272624529d-4, 2.89812192181749975d-4, &
      -1.44906096090874960d-4, 2.94615351118807751d-3, &
      -4.48891257791280696d-1, 1.87035260703064810d-1, 1.24690173802143589d-1, &
      -6.23450869010717876d-2, 3.70443771408077230d-2, &
      0.00000000000000000d+0, 0.00000000000000000d+0, 0.00000000000000000d+0, &
      0.00000000000000000d+0, 0.00000000000000000d+0, &
      0.00000000000000000d+0, -1.69406589450860068d-21, -5.08219768352580203d-21, &
      -8.47032947254300339d-22, -4.33680868994201774d-19, &
      -5.52498366688519482d-1, -1.88150728566460979d-1, -1.25433819044407258d-1, &
      6.27169095222036288d-2, 1.23798574626422572d-2, &
      -6.51783738721830708d-1, -4.26909843934396814d-3, -2.84606562622930718d-3, &
      1.42303281311465381d-3, 1.49521341773236269d-2, &
      -8.68012840637807648d-2, 1.06804545107300869d-3, 7.12030300715337720d-4, &
      -3.56015150357668914d-4, 2.77256453950402125d-3, &
      -1.73648344394211451d-1, 1.94003901717258427d-3, 1.29335934478171793d-3, &
      -6.46679672390859076d-4, 6.95957034219114347d-3, &
      2.80310059515267440d-1, -2.47572667013062691d-3, -1.65048444675374708d-3, &
      8.25242223376873541d-4, -2.82461774595950815d-2], [5,4,3])
    real(8) :: rr(1,4),gg(1,3,4),tt(1,4),ww(1),rd(1,4),gd(1,3,4),td(1,4),wd(1)
    real(8) :: fields(5,4),observed(5,4),direction(3),ev,err,limit
    integer :: case_index,component
    call quadrature(17)
    ww=1.d0
    do case_index=1,3
      direction=[0.d0,0.d0,1.d0]
      if(case_index==3) direction=[1.d0,2.d0,-3.d0]/sqrt(14.d0)
      fields(:,1)=[0.3d0,0.03d0,0.02d0,-0.01d0,0.2d0]
      do component=2,4
        fields(:,component)=polarization(case_index)*direction(component-1)*fields(:,1)
      end do
      rr(1,:)=fields(1,:)
      gg(1,:,:)=fields(2:4,:)
      tt(1,:)=fields(5,:)
      call SCAN$PAULI_BLOCK(1,nq,rr,gg,tt,ww,axes,weights,ev,rd,gd,td,wd)
      observed(1,:)=rd(1,:)
      observed(2:4,:)=gd(1,:,:)
      observed(5,:)=td(1,:)
      err=maxval(abs(observed-reference(:,:,case_index)))
      limit=1.d-10
      if(case_index==2) limit=1.d-7
      write(*,'(A,I2,A,ES18.10)') 'SCAN INDEPENDENT MC RESPONSE CASE ',case_index,' ERROR ',err
      if(abs(ev-eref(case_index))>1.d-12.or.err>limit) error stop 'Independent strong-spin SCAN reference'
    end do
  end subroutine
  subroutine vacuum_response
    real(8) :: base_r(np,4),base_g(np,3,4),base_t(np,4),factor,numeric,exact,relative
    integer :: scale
    call quadrature(3)
    do scale=2,8,2
      factor=10.d0**(-scale)
      base_r=factor*r0
      base_g=factor*g0
      base_t=factor*t0
      r=base_r
      g=base_g
      t=base_t
      ar=0.03d0*base_r
      ag=-0.04d0*base_g
      at=0.07d0*base_t
      call evaluate
      exact=sum(dr*ar)+sum(dg*ag)+sum(dt*at)
      h=1.d-4
      r=base_r+h*ar
      g=base_g+h*ag
      t=base_t+h*at
      call evaluate
      ep=e
      r=base_r-h*ar
      g=base_g-h*ag
      t=base_t-h*at
      call evaluate
      em=e
      numeric=(ep-em)/(2.d0*h)
      relative=abs(numeric-exact)/max(abs(numeric),abs(exact),1.d-25)
      write(*,'(A,I3,A,ES18.10)') 'SCAN DILUTE PAULI RESPONSE SCALE ',scale,' RELATIVE ERROR ',relative
      if(relative>2.d-6) error stop 'Dilute multicollinear SCAN response derivative'
    end do
  end subroutine

  subroutine quadrature(exactness)
    integer,intent(in) :: exactness
    if(allocated(axes)) deallocate(axes,weights)
    call paw_lebedev_size(exactness,nq)
    allocate(axes(3,nq),weights(nq))
    call paw_lebedev_get(exactness,nq,axes,weights)
  end subroutine
  subroutine perturb(h)
    real(8),intent(in) :: h
    r=r0+h*ar
    g=g0+h*ag
    t=t0+h*at
    x=x0
    w=w0
    c=c0
    aw=aw0
    if(k==4) then
      x=(1.d0+0.03d0*h)*x0
      w=(1.d0+0.03d0*h)*w0
      c=(1.d0+0.03d0*h)*c0
      aw=(1.d0+0.03d0*h)*aw0
    end if
  end subroutine
  subroutine evaluate
    call SKALA$EVALUATE_PAULI(np,na,nq,r,g,t,x,w,c,aw,sizes,axes,weights,e,dr,dg,dt,dx,dw,dc,daw)
    if(.not.ieee_is_finite(e).or.any(.not.ieee_is_finite(dr)).or.any(.not.ieee_is_finite(dg)) &
      .or.any(.not.ieee_is_finite(dt)).or.any(.not.ieee_is_finite(dx)).or.any(.not.ieee_is_finite(dw)) &
      .or.any(.not.ieee_is_finite(dc)).or.any(.not.ieee_is_finite(daw))) error stop 'Nonfinite Pauli XC'
  end subroutine
end program xc_spin_noncollinear
