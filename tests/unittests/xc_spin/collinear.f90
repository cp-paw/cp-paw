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
  integer :: i, ia, j, k
  logical :: available
  character(16) :: functional
  character(16) :: precision
  character(2048) :: model

  if(command_argument_count()<1) error stop 'Usage: xc_spin.x SCAN|SKALA [MODEL]'
  call get_command_argument(1,functional)
  if(functional/='SCAN'.and.functional/='SKALA') error stop 'Unknown XC functional'
  model=''
  if(command_argument_count()>=2) call get_command_argument(2,model)
  if(functional=='SKALA'.and.len_trim(model)==0) error stop 'Skala requires a model'
  precision='published'
  if(command_argument_count()>=3) call get_command_argument(3,precision)
  if(precision/='published'.and.precision/='diagnostic64') error stop 'Unknown model precision'
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
  write(*,'(A,A,A,ES18.10)') trim(functional),' CLOSED-SHELL SPIN SYMMETRY', &
                           ' SCALED ERROR ',err/scale
  if(err/scale>1.d-9) error stop 'Closed-shell spin-channel adjoints differ'

  ! Independent charge and spin directions exercise the two-channel chain rule.
  tolerance=5.d-8
  steps=[1.d-4,3.d-5]
  if(functional=='SKALA'.and.precision=='published') then
    ! The released export rounds intermediate features to Float32. A smaller
    ! displacement need not improve finite differences of that rounded energy.
    steps=[1.d-3,3.d-4]
    tolerance=2.d-5
  end if
  do k=1,2
    ar=0.07d0*r0
    ag=-0.04d0*g0
    at=0.09d0*t0
    if(k==2) then
      ar(:,2)=-ar(:,2)
      ag(:,:,2)=-ag(:,:,2)
      at(:,2)=-at(:,2)
    end if
    analytic=sum(dr0*ar)+sum(dg0*ag)+sum(dt0*at)
    if(k==2.and.abs(analytic)>1.d-9*scale) error stop 'Closed-shell magnetic derivative'
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
      write(*,'(A,A,I2,A,ES12.4,A,ES18.10)') trim(functional), &
        ' CLOSED-SHELL DIRECTION ',k,' STEP ',h,' SCALED FD ERROR ',err
      if(err>tolerance) error stop 'XC primitive directional derivative mismatch'
    end do
  end do
  call SKALA$FINALIZE
  call SKALA$SETL4('ON',.false.)
  write(*,'(A,A)') trim(functional),' CLOSED-SHELL SPIN CHECKS PASSED'
  call MPE$EXIT

contains

  subroutine evaluate
    call SKALA$EVALUATE(np,na,rho,g,tau,x,w,c,aw,sizes,e,dr,dg,dt,dx,dw,dc,daw)
    if(.not.ieee_is_finite(e).or.any(.not.ieee_is_finite(dr)).or. &
       any(.not.ieee_is_finite(dg)).or.any(.not.ieee_is_finite(dt)).or. &
       any(.not.ieee_is_finite(dx)).or.any(.not.ieee_is_finite(dw)).or. &
       any(.not.ieee_is_finite(dc)).or.any(.not.ieee_is_finite(daw))) &
      error stop 'Nonfinite XC energy or adjoint'
  end subroutine evaluate

end program xc_spin_collinear
