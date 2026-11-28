! SPDX-License-Identifier: GPL-3.0-or-later
program reject_spin_orbit
  use setup_module, only: this
  implicit none
  complex(8) :: dm(1,1,4,1)
  character(16) :: functional,reltype
  call get_command_argument(1,functional)
  if(functional/='SCAN'.and.functional/='SKALA') error stop 'Expected a native XC backend'
  call MPE$INIT
  call SKALA$SETCH('FUNCTIONAL',trim(functional))
  call SKALA$SETL4('MULTICOLLINEAR',.true.)
  ! Only the species metadata exist. No SOC radial or electronic calculation is run.
  call SETUP$SELECT('test')
  this%setting%trel=.true.
  this%setting%so=.true.
  call SETUP$GETCH('RELTYPE',reltype)
  if(reltype/='SPINORBIT') error stop 'Unexpected test setup type'
  call SETUP$UNSELECT()
  dm=(0.d0,0.d0)
  call SKALA$PREPARE(1,[1],1,4,dm)
  error stop 'Spin-orbit setup passed the guard'
end program
