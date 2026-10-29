program xc_spin_collinear
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: na=2, np=16
  real(8) :: rho(np,2), g(np,3,2), tau(np,2), x(np,3), w(np), c(na,3), aw(np)
  real(8) :: dr(np,2), dg(np,3,2), dt(np,2), dx(np,3), dw(np), dc(na,3), daw(np)
  real(8) :: r0(np,2), g0(np,3,2), t0(np,2), dr0(np,2), dg0(np,3,2), dt0(np,2)
  real(8) :: ar(np,2), ag(np,3,2), at(np,2), e, e0, ep, em, analytic, numeric
  real(8) :: angle, radius, err, scale, h, tolerance
  real(8) :: steps(2)
  integer(8) :: sizes(na)
  integer :: i, ia, j, k, ndirection
  logical :: available, open_shell
  character(16) :: functional
  character(16) :: precision
  character(16) :: spin_case
  character(40) :: label
  character(2048) :: model

  if(command_argument_count()<1.or.command_argument_count()>4) then
    write(*,'(A)') 'Usage: xc_spin.x SCAN|SKALA [MODEL [published|diagnostic64 [closed|open]]]'
    error stop 'Invalid XC spin test arguments'
  end if
  call get_command_argument(1,functional)
  if(functional/='SCAN'.and.functional/='SKALA') error stop 'Unknown XC functional'
  model=''
  if(command_argument_count()>=2) call get_command_argument(2,model)
  if(functional=='SKALA'.and.len_trim(model)==0) error stop 'Skala requires a model'
  precision='published'
  if(command_argument_count()>=3) call get_command_argument(3,precision)
  if(precision/='published'.and.precision/='diagnostic64') error stop 'Unknown model precision'
  spin_case='closed'
  if(command_argument_count()>=4) call get_command_argument(4,spin_case)
  if(spin_case/='closed'.and.spin_case/='open') error stop 'Unknown spin case'
  open_shell=spin_case=='open'
  label=trim(functional)//' CLOSED-SHELL'
  if(open_shell) label=trim(functional)//' OPEN-SHELL'
  call MPE$INIT
  call SKALA$SETCH('FUNCTIONAL',trim(functional))
  call SKALA$SETL4('ON',.true.)
  call SKALA$SETL4('APPLY',.true.)
  call SKALA$SETL4('CHECK',.false.)
  if(functional=='SKALA') then
    call SKALA$SETCH('MODEL',trim(model))
    call SKALA$SETCH('DEVICE','CPU')
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
    radius=0.35d0+0.03d0*real(mod(j,3),8)
    x(i,:)=c(ia,:)+[radius*cos(angle),radius*sin(angle),0.02d0*real(j-4,8)]
    rho(i,:)=0.25d0*exp(-radius)
    g(i,:,1)=[-0.12d0*cos(angle),-0.12d0*sin(angle),0.01d0*sin(angle)]
    g(i,:,2)=g(i,:,1)
    tau(i,:)=0.18d0*rho(i,:)
    w(i)=0.04d0+0.001d0*real(j,8)
    aw(i)=0.05d0
  end do
  if(open_shell) then
    rho(:,1)=1.35d0*rho(:,1)
    rho(:,2)=0.65d0*rho(:,2)
    g(:,:,1)=1.25d0*g(:,:,1)
    g(:,:,2)=0.75d0*g(:,:,2)
    tau(:,1)=1.15d0*tau(:,1)
    tau(:,2)=0.85d0*tau(:,2)
  end if
  r0=rho
  g0=g
  t0=tau
  call evaluate
  e0=e
  dr0=dr
  dg0=dg
  dt0=dt
  if(abs(e0)<1.d-12.or.maxval(abs(dr0))<1.d-12) error stop 'Degenerate XC result'
  err=max(maxval(abs(dr(:,1)-dr(:,2))),maxval(abs(dg(:,:,1)-dg(:,:,2))), &
          maxval(abs(dt(:,1)-dt(:,2))))
  scale=max(1.d0,maxval(abs(dr)),maxval(abs(dg)),maxval(abs(dt)))
  if(open_shell) then
    if(err/scale<1.d-9) error stop 'Open-shell XC adjoints lost their spin dependence'
  else
    write(*,'(A,A,ES18.10)') trim(label),' SPIN SYMMETRY SCALED ERROR ',err/scale
    if(err/scale>1.d-9) error stop 'Closed-shell spin-channel adjoints differ'
  end if
  call spin_swap(trim(label))

  ! Independent charge and spin directions exercise the two-channel chain rule.
  tolerance=5.d-8
  steps=[1.d-4,3.d-5]
  if(functional=='SKALA'.and.precision=='published') then
    ! The released export rounds intermediate features to Float32. A smaller
    ! displacement need not improve finite differences of that rounded energy.
    steps=[1.d-3,3.d-4]
    tolerance=2.d-5
  end if
  ndirection=2
  if(open_shell) ndirection=6
  do k=1,ndirection
    ar=0.07d0*r0
    ag=-0.04d0*g0
    at=0.09d0*t0
    if(k==2) then
      ar(:,2)=-ar(:,2)
      ag(:,:,2)=-ag(:,:,2)
      at(:,2)=-at(:,2)
    end if
    if(k>2) then
      ar=0.d0
      ag=0.d0
      at=0.d0
      select case(k)
      case(3)
        ar(:,1)=0.07d0*r0(:,1)
      case(4)
        ag(:,:,2)=-0.04d0*g0(:,:,2)
      case(5)
        at(:,1)=0.09d0*t0(:,1)
      case(6)
        at(:,2)=0.11d0*t0(:,2)
      end select
    end if
    analytic=sum(dr0*ar)+sum(dg0*ag)+sum(dt0*at)
    if(.not.open_shell.and.k==2.and.abs(analytic)>1.d-9*scale) &
      error stop 'Closed-shell magnetic derivative'
    do j=1,size(steps)
      h=steps(j)
      rho=r0+h*ar
      g=g0+h*ag
      tau=t0+h*at
      call evaluate
      ep=e
      rho=r0-h*ar
      g=g0-h*ag
      tau=t0-h*at
      call evaluate
      em=e
      numeric=(ep-em)/(2.d0*h)
      err=abs(numeric-analytic)/max(1.d0,abs(analytic))
      write(*,'(A,A,I2,A,ES12.4,A,ES18.10)') trim(label), &
        ' DIRECTION ',k,' STEP ',h,' SCALED FD ERROR ',err
      if(err>tolerance) error stop 'XC primitive directional derivative mismatch'
    end do
  end do
  if(open_shell) then
    rho=0.d0
    g=0.d0
    tau=0.d0
    rho(:,1)=sum(r0,dim=2)
    g(:,:,1)=sum(g0,dim=3)
    tau(:,1)=sum(t0,dim=2)
    call evaluate
    call spin_swap(trim(functional)//' EMPTY SPIN CHANNEL')
  end if
  call SKALA$FINALIZE
  call SKALA$SETL4('ON',.false.)
  write(*,'(A,A)') trim(label),' SPIN CHECKS PASSED'
  call MPE$EXIT

contains

  subroutine spin_swap(description)
    character(*), intent(in) :: description
    real(8) :: sr(np,2), sg(np,3,2), st(np,2), se
    real(8) :: sdr(np,2), sdg(np,3,2), sdt(np,2)
    real(8) :: sdx(np,3), sdw(np), sdc(na,3), sdaw(np), difference, norm, bound
    sr=rho
    sg=g
    st=tau
    se=e
    sdr=dr
    sdg=dg
    sdt=dt
    sdx=dx
    sdw=dw
    sdc=dc
    sdaw=daw
    rho=sr(:,2:1:-1)
    g=sg(:,:,2:1:-1)
    tau=st(:,2:1:-1)
    call evaluate
    difference=max(abs(e-se),maxval(abs(dr-sdr(:,2:1:-1))), &
                   maxval(abs(dg-sdg(:,:,2:1:-1))),maxval(abs(dt-sdt(:,2:1:-1))), &
                   maxval(abs(dx-sdx)),maxval(abs(dw-sdw)),maxval(abs(dc-sdc)), &
                   maxval(abs(daw-sdaw)))
    norm=max(1.d0,abs(se),maxval(abs(sdr)),maxval(abs(sdg)),maxval(abs(sdt)), &
             maxval(abs(sdx)),maxval(abs(sdw)),maxval(abs(sdc)),maxval(abs(sdaw)))
    bound=5.d-12
    if(functional=='SKALA') bound=5.d-10
    if(functional=='SKALA'.and.precision=='published') bound=2.d-7
    write(*,'(A,A,ES18.10)') description,' SPIN SWAP SCALED ERROR ',difference/norm
    if(difference/norm>bound) error stop 'XC energy or adjoint violates spin interchange'
    rho=sr
    g=sg
    tau=st
    e=se
    dr=sdr
    dg=sdg
    dt=sdt
    dx=sdx
    dw=sdw
    dc=sdc
    daw=sdaw
  end subroutine spin_swap

  subroutine evaluate
    call SKALA$EVALUATE(np,na,rho,g,tau,x,w,c,aw,sizes,e,dr,dg,dt,dx,dw,dc,daw)
    if(.not.ieee_is_finite(e).or.any(.not.ieee_is_finite(dr)).or. &
       any(.not.ieee_is_finite(dg)).or.any(.not.ieee_is_finite(dt)).or. &
       any(.not.ieee_is_finite(dx)).or.any(.not.ieee_is_finite(dw)).or. &
       any(.not.ieee_is_finite(dc)).or.any(.not.ieee_is_finite(daw))) &
      error stop 'Nonfinite XC energy or adjoint'
  end subroutine evaluate

end program xc_spin_collinear
