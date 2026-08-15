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
  real(real64), target :: grid_coord_deriv(npoint, 3), grid_weight_deriv(npoint)
  real(real64), target :: atom_coord_deriv(natom, 3), atomic_weight_deriv(npoint)
  integer(int64), target :: atomic_grid_sizes(natom)
  real(real64), parameter :: fd_step = 1.0e-4_real64, fd_tolerance = 2.0e-5_real64
  real(real64), parameter :: repeat_tolerance = 2.0e-8_real64
  real(real64) :: energy, angle, radius, eplus, eminus, original, max_error
  real(real64) :: repeat_energy, repeat_error, repeat_values(7)
  real(real64) :: translation_error, translation_gradient(3)
  real(real64) :: analytic(7)
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
    & grad_deriv, kin_deriv, status, message, grid_coord_deriv, grid_weight_deriv, &
    & atom_coord_deriv, atomic_weight_deriv)
  if (status /= 0) then
    write (*, '(a)') "Skala evaluation failed: " // trim(message)
    stop 4
  end if
  if (.not. ieee_is_finite(energy) .or. .not. all(ieee_is_finite(density_deriv)) &
      .or. .not. all(ieee_is_finite(grad_deriv)) .or. .not. all(ieee_is_finite(kin_deriv))) then
    write (*, '(a)') "Skala evaluation returned a non-finite value"
    stop 5
  end if
  if (.not. all(ieee_is_finite(grid_coord_deriv)) &
      .or. .not. all(ieee_is_finite(grid_weight_deriv)) &
      .or. .not. all(ieee_is_finite(atom_coord_deriv)) &
      .or. .not. all(ieee_is_finite(atomic_weight_deriv))) then
    write (*, '(a)') "Skala coordinate or weight gradient is non-finite"
    stop 5
  end if

  analytic = [density_deriv(3, 2), grad_deriv(4, 2, 1), kin_deriv(5, 2), &
    & grid_coord_deriv(6, 3), grid_weight_deriv(7), atom_coord_deriv(2, 3), &
    & atomic_weight_deriv(9)]
  repeat_energy = energy
  call paw_skala_evaluate(model, density, grad, kin, grid_coords, grid_weights, &
    & atom_coords, atomic_grid_weights, atomic_grid_sizes, energy, density_deriv, &
    & grad_deriv, kin_deriv, status, message, grid_coord_deriv, grid_weight_deriv, &
    & atom_coord_deriv, atomic_weight_deriv)
  if (status /= 0) then
    write (*, '(a)') "Skala repeat evaluation failed: " // trim(message)
    stop 4
  end if
  repeat_values = [density_deriv(3, 2), grad_deriv(4, 2, 1), kin_deriv(5, 2), &
    & grid_coord_deriv(6, 3), grid_weight_deriv(7), atom_coord_deriv(2, 3), &
    & atomic_weight_deriv(9)]
  repeat_error = max(abs(energy - repeat_energy), maxval(abs(repeat_values - analytic)))
  write (*, '(a,es12.4)') "SKALA_BRIDGE_REPEAT max absolute difference=", repeat_error
  if (repeat_error > repeat_tolerance) then
    write (*, '(a)') "Skala repeated evaluation is not reproducible"
    stop 8
  end if
  max_error = 0.0_real64
  translation_gradient = sum(grid_coord_deriv, dim=1) + sum(atom_coord_deriv, dim=1)
  translation_error = maxval(abs(translation_gradient))
  write (*, '(a,3es16.8,a,es12.4)') "SKALA_BRIDGE_TRANSLATION gradient=", &
    & translation_gradient, " max_error=", translation_error
  if (translation_error > 1.0e-10_real64) then
    write (*, '(a)') "Skala coordinate gradients violate translation invariance"
    stop 7
  end if

  original = density(3, 2)
  density(3, 2) = original + fd_step
  call evaluate_energy(eplus)
  density(3, 2) = original - fd_step
  call evaluate_energy(eminus)
  density(3, 2) = original
  call compare_gradient("density", analytic(1), eplus, eminus)

  original = grad(4, 2, 1)
  grad(4, 2, 1) = original + fd_step
  call evaluate_energy(eplus)
  grad(4, 2, 1) = original - fd_step
  call evaluate_energy(eminus)
  grad(4, 2, 1) = original
  call compare_gradient("grad", analytic(2), eplus, eminus)

  original = kin(5, 2)
  kin(5, 2) = original + fd_step
  call evaluate_energy(eplus)
  kin(5, 2) = original - fd_step
  call evaluate_energy(eminus)
  kin(5, 2) = original
  call compare_gradient("kin", analytic(3), eplus, eminus)

  original = grid_coords(6, 3)
  grid_coords(6, 3) = original + fd_step
  call evaluate_energy(eplus)
  grid_coords(6, 3) = original - fd_step
  call evaluate_energy(eminus)
  grid_coords(6, 3) = original
  call compare_gradient("grid_coords", analytic(4), eplus, eminus)

  original = grid_weights(7)
  grid_weights(7) = original + fd_step
  call evaluate_energy(eplus)
  grid_weights(7) = original - fd_step
  call evaluate_energy(eminus)
  grid_weights(7) = original
  call compare_gradient("grid_weights", analytic(5), eplus, eminus)

  original = atom_coords(2, 3)
  atom_coords(2, 3) = original + fd_step
  call evaluate_energy(eplus)
  atom_coords(2, 3) = original - fd_step
  call evaluate_energy(eminus)
  atom_coords(2, 3) = original
  call compare_gradient("atom_coords", analytic(6), eplus, eminus)

  original = atomic_grid_weights(9)
  atomic_grid_weights(9) = original + fd_step
  call evaluate_energy(eplus)
  atomic_grid_weights(9) = original - fd_step
  call evaluate_energy(eminus)
  atomic_grid_weights(9) = original
  call compare_gradient("atomic_grid_weights", analytic(7), eplus, eminus)

  write (*, '(a,es24.16)') "SKALA_BRIDGE_SMOKE energy=", energy
  write (*, '(a,es12.4)') "SKALA_BRIDGE_SMOKE max gradient error=", max_error
  write (*, '(a)') "SKALA_BRIDGE_SMOKE passed"

contains

  subroutine evaluate_energy(current_energy)
    real(real64), intent(out) :: current_energy
    call paw_skala_evaluate(model, density, grad, kin, grid_coords, grid_weights, &
      & atom_coords, atomic_grid_weights, atomic_grid_sizes, current_energy, density_deriv, &
      & grad_deriv, kin_deriv, status, message)
    if (status /= 0) then
      write (*, '(a)') "Skala finite-difference evaluation failed: " // trim(message)
      stop 6
    end if
  end subroutine evaluate_energy

  subroutine compare_gradient(label, expected, plus_energy, minus_energy)
    character(len=*), intent(in) :: label
    real(real64), intent(in) :: expected, plus_energy, minus_energy
    real(real64) :: error, finite_difference, scale
    finite_difference = (plus_energy - minus_energy) / (2.0_real64 * fd_step)
    error = abs(expected - finite_difference)
    scale = max(1.0_real64, abs(expected), abs(finite_difference))
    max_error = max(max_error, error / scale)
    write (*, '(a,1x,a,2(a,es16.8),a,es12.4)') "SKALA_BRIDGE_FD", trim(label), &
      & " analytic=", expected, " finite_difference=", finite_difference, &
      & " relative_error=", error / scale
    if (error > fd_tolerance * scale) then
      write (*, '(a)') "Skala feature-gradient finite-difference check failed"
      stop 7
    end if
  end subroutine compare_gradient
end program skala_bridge_smoke
