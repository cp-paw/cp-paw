! SPDX-License-Identifier: GPL-3.0-or-later

program skala_bridge_smoke
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use paw_skala_bridge, only : paw_skala_evaluate, paw_skala_load, paw_skala_model
  implicit none

  integer, parameter :: natom = 2, points_per_atom = 8, npoint = natom * points_per_atom
  type(paw_skala_model) :: model
  real(real64), target :: density(npoint, 2), grad(npoint, 3, 2), kin(npoint, 2)
  real(real64), target :: grid_coords(npoint, 3), grid_weights(npoint)
  real(real64), target :: atom_coords(natom, 3), atomic_grid_weights(npoint)
  real(real64), target :: density_deriv(npoint, 2), grad_deriv(npoint, 3, 2)
  real(real64), target :: kin_deriv(npoint, 2)
  integer(int64), target :: atomic_grid_sizes(natom)
  real(real64) :: energy, angle, radius
  integer :: device_index, i, iatom, ipoint, status
  character(len=1024) :: model_path, device, device_index_text, message

  if (command_argument_count() < 1) then
    write (*, '(a)') "usage: cppaw_skala_smoke MODEL [cpu|cuda] [DEVICE_INDEX]"
    stop 2
  end if
  call get_command_argument(1, model_path)
  device = "cpu"
  device_index = 0
  if (command_argument_count() >= 2) call get_command_argument(2, device)
  if (command_argument_count() >= 3) then
    call get_command_argument(3, device_index_text)
    read (device_index_text, *, iostat=status) device_index
    if (status /= 0) then
      write (*, '(a)') "invalid CUDA device index"
      stop 2
    end if
  end if

  atom_coords = 0.0_real64
  atom_coords(2, 3) = 1.4_real64
  atomic_grid_sizes = points_per_atom
  do iatom = 1, natom
    do i = 1, points_per_atom
      ipoint = (iatom - 1) * points_per_atom + i
      angle = 2.0_real64 * acos(-1.0_real64) * real(i - 1, real64) / &
        & real(points_per_atom, real64)
      radius = 0.35_real64 + 0.03_real64 * real(mod(i, 3), real64)
      grid_coords(ipoint, :) = atom_coords(iatom, :)
      grid_coords(ipoint, 1) = grid_coords(ipoint, 1) + radius * cos(angle)
      grid_coords(ipoint, 2) = grid_coords(ipoint, 2) + radius * sin(angle)
      grid_coords(ipoint, 3) = grid_coords(ipoint, 3) + 0.02_real64 * real(i - 4, real64)
      density(ipoint, 1) = 0.25_real64 * exp(-radius)
      density(ipoint, 2) = density(ipoint, 1)
      grad(ipoint, :, :) = 0.0_real64
      grad(ipoint, 1, :) = -0.12_real64 * cos(angle)
      grad(ipoint, 2, :) = -0.12_real64 * sin(angle)
      kin(ipoint, :) = 0.18_real64 * density(ipoint, :)
      grid_weights(ipoint) = 0.04_real64
      atomic_grid_weights(ipoint) = 0.05_real64
    end do
  end do

  call paw_skala_load(model, trim(model_path), trim(device), device_index, status, message)
  if (status /= 0) then
    write (*, '(a)') "Skala load failed: " // trim(message)
    stop 3
  end if
  call paw_skala_evaluate(model, density, grad, kin, grid_coords, grid_weights, &
    & atom_coords, atomic_grid_weights, atomic_grid_sizes, energy, density_deriv, &
    & grad_deriv, kin_deriv, status, message)
  if (status /= 0) then
    write (*, '(a)') "Skala evaluation failed: " // trim(message)
    stop 4
  end if
  if (.not. ieee_is_finite(energy) .or. .not. all(ieee_is_finite(density_deriv)) &
      .or. .not. all(ieee_is_finite(grad_deriv)) .or. .not. all(ieee_is_finite(kin_deriv))) then
    write (*, '(a)') "Skala evaluation returned a non-finite value"
    stop 5
  end if
  write (*, '(a,es24.16)') "SKALA_BRIDGE_SMOKE energy=", energy
  write (*, '(a)') "SKALA_BRIDGE_SMOKE passed"
end program skala_bridge_smoke
