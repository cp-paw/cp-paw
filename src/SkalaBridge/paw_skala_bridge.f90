! SPDX-License-Identifier: GPL-3.0-or-later

module paw_skala_bridge
  use, intrinsic :: iso_c_binding, only : c_associated, c_char, c_double, c_int, &
    & c_int64_t, c_null_char, c_null_ptr, c_ptr
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use ftorch, only : ftorch_int, torch_delete, torch_kCPU, torch_kCUDA, &
    & torch_tensor, torch_tensor_from_array
  implicit none
  private

  integer, parameter :: error_capacity = 2048
  integer, parameter :: nfeatures = 9

  type, public :: paw_skala_model
    private
    type(c_ptr) :: handle = c_null_ptr
    integer(c_int) :: device_type = torch_kCPU
    integer(c_int) :: device_index = -1
    logical :: features(nfeatures) = .false.
  contains
    final :: finalize_model
  end type paw_skala_model

  public :: paw_skala_load, paw_skala_evaluate, paw_skala_release
  public :: paw_skala_cuda_available
  public :: paw_skala_device_cpu, paw_skala_device_cuda

contains

  logical function paw_skala_cuda_available()
    interface
      function cuda_available_c() result(available) &
          bind(c, name="cppaw_skala_cuda_available")
        import :: c_int
        integer(c_int) :: available
      end function cuda_available_c
    end interface
    paw_skala_cuda_available = cuda_available_c() /= 0
  end function paw_skala_cuda_available

  integer function paw_skala_device_cpu()
    paw_skala_device_cpu = torch_kCPU
  end function paw_skala_device_cpu

  integer function paw_skala_device_cuda()
    paw_skala_device_cuda = torch_kCUDA
  end function paw_skala_device_cuda

  subroutine paw_skala_load(model, path, device, device_index, status, message)
    type(paw_skala_model), intent(inout) :: model
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: device
    integer, intent(in), optional :: device_index
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    interface
      function model_load_c(path_c, device_type, index, features, error, capacity) &
          result(handle) bind(c, name="cppaw_skala_model_load")
        import :: c_char, c_int, c_ptr
        character(kind=c_char), intent(in) :: path_c(*)
        integer(c_int), value :: device_type, index, capacity
        integer(c_int), intent(out) :: features(*)
        character(kind=c_char), intent(out) :: error(*)
        type(c_ptr) :: handle
      end function model_load_c
    end interface

    character(kind=c_char) :: error(error_capacity)
    integer(c_int) :: feature_flags(nfeatures)

    call paw_skala_release(model)
    select case (trim(lowercase(device)))
    case ("cpu")
      model%device_type = torch_kCPU
      model%device_index = -1
    case ("cuda")
      model%device_type = torch_kCUDA
      model%device_index = 0
      if (present(device_index)) model%device_index = max(0, device_index)
    case default
      status = 1
      message = "Skala device must be 'cpu' or 'cuda'"
      return
    end select

    error = c_null_char
    feature_flags = 0
    model%handle = model_load_c(c_string(path), model%device_type, &
      & model%device_index, feature_flags, error, error_capacity)
    if (.not. c_associated(model%handle)) then
      status = 1
      call copy_c_message(error, message)
      return
    end if
    model%features = feature_flags /= 0
    if (.not. all(model%features)) then
      status = 1
      message = "Skala model does not provide the complete protocol-v2 feature set"
      call paw_skala_release(model)
      return
    end if
    status = 0
    message = ""
  end subroutine paw_skala_load

  subroutine paw_skala_evaluate(model, density, grad, kin, grid_coords, grid_weights, &
      & atom_coords, atomic_grid_weights, atomic_grid_sizes, energy, density_deriv, &
      & grad_deriv, kin_deriv, status, message, grid_coord_deriv, grid_weight_deriv, &
      & atom_coord_deriv, atomic_weight_deriv)
    type(paw_skala_model), intent(inout) :: model
    real(real64), intent(in) :: density(:,:), grad(:,:,:), kin(:,:)
    real(real64), intent(in), target :: grid_coords(:,:), grid_weights(:)
    real(real64), intent(in), target :: atom_coords(:,:), atomic_grid_weights(:)
    integer(int64), intent(in), target :: atomic_grid_sizes(:)
    real(real64), intent(out) :: energy
    real(real64), intent(out), target :: density_deriv(:,:), grad_deriv(:,:,:), kin_deriv(:,:)
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(real64), intent(out), optional, target :: grid_coord_deriv(:,:), grid_weight_deriv(:)
    real(real64), intent(out), optional, target :: atom_coord_deriv(:,:), atomic_weight_deriv(:)

    interface
      function dict_create_c() result(handle) bind(c, name="cppaw_skala_dict_create")
        import :: c_ptr
        type(c_ptr) :: handle
      end function dict_create_c
      subroutine dict_release_c(handle) bind(c, name="cppaw_skala_dict_release")
        import :: c_ptr
        type(c_ptr), value :: handle
      end subroutine dict_release_c
      function dict_insert_c(handle, key, tensor, error, capacity) result(rc) &
          bind(c, name="cppaw_skala_dict_insert")
        import :: c_char, c_int, c_ptr
        type(c_ptr), value :: handle, tensor
        character(kind=c_char), intent(in) :: key(*)
        character(kind=c_char), intent(out) :: error(*)
        integer(c_int), value :: capacity
        integer(c_int) :: rc
      end function dict_insert_c
      function evaluate_c(model_handle, dict_handle, energy_c, d_density, d_grad, &
          d_kin, d_grid_coords, d_grid_weights, d_atom_coords, d_atomic_weights, &
          error, capacity) result(rc) bind(c, name="cppaw_skala_model_evaluate")
        import :: c_char, c_double, c_int, c_ptr
        type(c_ptr), value :: model_handle, dict_handle
        real(c_double), intent(out) :: energy_c
        type(c_ptr), intent(out) :: d_density, d_grad, d_kin, d_grid_coords
        type(c_ptr), intent(out) :: d_grid_weights, d_atom_coords, d_atomic_weights
        character(kind=c_char), intent(out) :: error(*)
        integer(c_int), value :: capacity
        integer(c_int) :: rc
      end function evaluate_c
      function tensor_copy_c(tensor, destination, count, error, capacity) result(rc) &
          bind(c, name="cppaw_skala_tensor_copy_double")
        import :: c_char, c_double, c_int, c_int64_t, c_ptr
        type(c_ptr), value :: tensor
        real(c_double), intent(out) :: destination(*)
        integer(c_int64_t), value :: count
        character(kind=c_char), intent(out) :: error(*)
        integer(c_int), value :: capacity
        integer(c_int) :: rc
      end function tensor_copy_c
      subroutine tensor_release_c(tensor) bind(c, name="cppaw_skala_tensor_release")
        import :: c_ptr
        type(c_ptr), value :: tensor
      end subroutine tensor_release_c
    end interface

    integer :: i, ipoint, npoint, natom, max_grid_size
    integer(c_int) :: rc
    integer(ftorch_int), parameter :: layout1(1) = [1]
    integer(ftorch_int), parameter :: layout2(2) = [1, 2]
    integer(ftorch_int), parameter :: layout3(3) = [1, 2, 3]
    real(real64), allocatable, target :: density_t_data(:,:), grad_t_data(:,:,:), kin_t_data(:,:)
    real(real64), allocatable, target :: grid_coord_t_data(:,:), atom_coord_t_data(:,:)
    integer(int64), allocatable, target :: bound_shape(:,:)
    type(torch_tensor) :: tensors(9)
    type(c_ptr) :: dict_handle
    type(c_ptr) :: derivative_handles(7)
    character(kind=c_char) :: error(error_capacity)

    energy = 0.0_real64
    density_deriv = 0.0_real64
    grad_deriv = 0.0_real64
    kin_deriv = 0.0_real64
    if (present(grid_coord_deriv)) grid_coord_deriv = 0.0_real64
    if (present(grid_weight_deriv)) grid_weight_deriv = 0.0_real64
    if (present(atom_coord_deriv)) atom_coord_deriv = 0.0_real64
    if (present(atomic_weight_deriv)) atomic_weight_deriv = 0.0_real64
    status = 1
    message = ""

    if (.not. c_associated(model%handle)) then
      message = "Skala model is not loaded"
      return
    end if
    npoint = size(density, 1)
    natom = size(atom_coords, 1)
    if (size(density, 2) /= 2 .or. size(kin, 1) /= npoint .or. size(kin, 2) /= 2 &
        .or. size(grad, 1) /= npoint .or. size(grad, 2) /= 3 .or. size(grad, 3) /= 2 &
        .or. size(grid_coords, 1) /= npoint .or. size(grid_coords, 2) /= 3 &
        .or. size(grid_weights) /= npoint .or. size(atomic_grid_weights) /= npoint &
        .or. size(atom_coords, 2) /= 3 .or. size(atomic_grid_sizes) /= natom &
        .or. sum(atomic_grid_sizes) /= npoint) then
      message = "inconsistent Skala feature dimensions"
      return
    end if
    if (any(atomic_grid_sizes <= 0_int64)) then
      message = "each Skala atom block must contain at least one grid point"
      return
    end if
    max_grid_size = int(maxval(atomic_grid_sizes))

    allocate(density_t_data(2, npoint), grad_t_data(2, 3, npoint), &
      & kin_t_data(2, npoint), grid_coord_t_data(3, npoint), &
      & atom_coord_t_data(3, natom), bound_shape(max_grid_size, 0))
    do ipoint = 1, npoint
      density_t_data(:, ipoint) = density(ipoint, :)
      kin_t_data(:, ipoint) = kin(ipoint, :)
      do i = 1, 3
        grad_t_data(:, i, ipoint) = grad(ipoint, i, :)
      end do
    end do

    call torch_tensor_from_array(tensors(1), density_t_data, layout2, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(2), grad_t_data, layout3, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(3), kin_t_data, layout2, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(4), grid_coords, layout2, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(5), grid_weights, layout1, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(6), atom_coords, layout2, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(7), atomic_grid_weights, layout1, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(8), atomic_grid_sizes, layout1, &
      & model%device_type, model%device_index)
    call torch_tensor_from_array(tensors(9), bound_shape, layout2, &
      & model%device_type, model%device_index)

    dict_handle = dict_create_c()
    call insert_feature("density", tensors(1))
    if (status == 0) call insert_feature("grad", tensors(2))
    if (status == 0) call insert_feature("kin", tensors(3))
    if (status == 0) call insert_feature("grid_coords", tensors(4))
    if (status == 0) call insert_feature("grid_weights", tensors(5))
    if (status == 0) call insert_feature("coarse_0_atomic_coords", tensors(6))
    if (status == 0) call insert_feature("atomic_grid_weights", tensors(7))
    if (status == 0) call insert_feature("atomic_grid_sizes", tensors(8))
    if (status == 0) call insert_feature("atomic_grid_size_bound_shape", tensors(9))
    if (status /= 0) then
      call cleanup()
      return
    end if

    derivative_handles = c_null_ptr
    error = c_null_char
    rc = evaluate_c(model%handle, dict_handle, energy, derivative_handles(1), &
      & derivative_handles(2), derivative_handles(3), derivative_handles(4), &
      & derivative_handles(5), derivative_handles(6), derivative_handles(7), &
      & error, error_capacity)
    if (rc /= 0) then
      call copy_c_message(error, message)
      call cleanup()
      return
    end if

    call copy_gradient(derivative_handles(1), density_deriv, int(size(density_deriv), int64))
    if (status == 0) call copy_gradient(derivative_handles(2), grad_deriv, int(size(grad_deriv), int64))
    if (status == 0) call copy_gradient(derivative_handles(3), kin_deriv, int(size(kin_deriv), int64))
    if (status == 0 .and. present(grid_coord_deriv)) then
      call copy_gradient(derivative_handles(4), grid_coord_t_data, int(size(grid_coord_t_data), int64))
      if (status == 0) grid_coord_deriv = transpose(grid_coord_t_data)
    end if
    if (status == 0 .and. present(grid_weight_deriv)) &
      & call copy_gradient(derivative_handles(5), grid_weight_deriv, int(size(grid_weight_deriv), int64))
    if (status == 0 .and. present(atom_coord_deriv)) then
      call copy_gradient(derivative_handles(6), atom_coord_t_data, int(size(atom_coord_t_data), int64))
      if (status == 0) atom_coord_deriv = transpose(atom_coord_t_data)
    end if
    if (status == 0 .and. present(atomic_weight_deriv)) &
      & call copy_gradient(derivative_handles(7), atomic_weight_deriv, int(size(atomic_weight_deriv), int64))
    call cleanup()

  contains

    subroutine insert_feature(key, tensor)
      character(len=*), intent(in) :: key
      type(torch_tensor), intent(in) :: tensor
      error = c_null_char
      rc = dict_insert_c(dict_handle, c_string(key), tensor%p, error, error_capacity)
      if (rc == 0) then
        status = 0
      else
        status = 1
        call copy_c_message(error, message)
      end if
    end subroutine insert_feature

    subroutine copy_gradient(handle, destination, count)
      type(c_ptr), intent(in) :: handle
      real(real64), intent(out) :: destination(*)
      integer(int64), intent(in) :: count
      error = c_null_char
      rc = tensor_copy_c(handle, destination, int(count, c_int64_t), error, error_capacity)
      if (rc /= 0) then
        status = 1
        call copy_c_message(error, message)
      else
        status = 0
      end if
    end subroutine copy_gradient

    subroutine cleanup()
      integer :: itensor
      if (c_associated(dict_handle)) call dict_release_c(dict_handle)
      dict_handle = c_null_ptr
      do itensor = 1, size(derivative_handles)
        if (c_associated(derivative_handles(itensor))) &
          & call tensor_release_c(derivative_handles(itensor))
      end do
      call torch_delete(tensors)
    end subroutine cleanup

  end subroutine paw_skala_evaluate

  subroutine paw_skala_release(model)
    type(paw_skala_model), intent(inout) :: model
    interface
      subroutine model_release_c(handle) bind(c, name="cppaw_skala_model_release")
        import :: c_ptr
        type(c_ptr), value :: handle
      end subroutine model_release_c
    end interface
    if (c_associated(model%handle)) call model_release_c(model%handle)
    model%handle = c_null_ptr
    model%features = .false.
  end subroutine paw_skala_release

  subroutine finalize_model(model)
    type(paw_skala_model), intent(inout) :: model
    call paw_skala_release(model)
  end subroutine finalize_model

  function c_string(value) result(converted)
    character(len=*), intent(in) :: value
    character(kind=c_char) :: converted(len_trim(value) + 1)
    integer :: i
    do i = 1, len_trim(value)
      converted(i) = value(i:i)
    end do
    converted(len_trim(value) + 1) = c_null_char
  end function c_string

  subroutine copy_c_message(source, destination)
    character(kind=c_char), intent(in) :: source(:)
    character(len=*), intent(out) :: destination
    integer :: i
    destination = ""
    do i = 1, min(size(source), len(destination))
      if (source(i) == c_null_char) exit
      destination(i:i) = source(i)
    end do
  end subroutine copy_c_message

  pure function lowercase(value) result(converted)
    character(len=*), intent(in) :: value
    character(len=len(value)) :: converted
    integer :: i, code
    converted = value
    do i = 1, len(value)
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) &
        & converted(i:i) = achar(code + iachar('a') - iachar('A'))
    end do
  end function lowercase

end module paw_skala_bridge
