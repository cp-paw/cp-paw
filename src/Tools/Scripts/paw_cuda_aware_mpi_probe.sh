#!/usr/bin/env bash
set -euo pipefail

RANKS=${RANKS:-2}
N=${N:-1024}
WORKDIR=${WORKDIR:-$(mktemp -d "${TMPDIR:-/tmp}/cppaw-cuda-aware-mpi.XXXXXX")}
KEEP_WORKDIR=${KEEP_WORKDIR:-no}
MPIFORT=${MPIFORT:-$(command -v mpifort 2>/dev/null || true)}
MPIRUN=${MPIRUN:-$(command -v mpirun 2>/dev/null || true)}
FCFLAGS=${FCFLAGS:-"-acc"}
MPI_ARGS=${MPI_ARGS:-}

cleanup() {
  case "${KEEP_WORKDIR}" in
    yes|true|1) ;;
    *) rm -rf "${WORKDIR}" ;;
  esac
}
trap cleanup EXIT

if [[ -z "${MPIFORT}" ]]; then
  echo "cuda_aware_mpi_probe=skip reason=no_mpifort"
  exit 0
fi
if [[ -z "${MPIRUN}" ]]; then
  echo "cuda_aware_mpi_probe=skip reason=no_mpirun"
  exit 0
fi

mkdir -p "${WORKDIR}"
src="${WORKDIR}/cuda_aware_mpi_probe.f90"
exe="${WORKDIR}/cuda_aware_mpi_probe.x"

cat > "${src}" <<'EOF'
program cuda_aware_mpi_probe
  use mpi
  implicit none
  integer :: ierr, rank, nranks, n, i
  real(8), allocatable :: send(:), recv(:)
  real(8) :: expected, local_err, global_err
  character(len=32) :: arg

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

  n = 1024
  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) read(arg, *) n
  allocate(send(n), recv(n))

  do i = 1, n
    send(i) = real(rank + 1, kind=8)
    recv(i) = -1.0d0
  end do

!$acc data copyin(send(1:n)) copyout(recv(1:n))
!$acc host_data use_device(send, recv)
  call MPI_Allreduce(send, recv, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!$acc end host_data
!$acc end data

  expected = 0.5d0 * real(nranks * (nranks + 1), kind=8)
  local_err = maxval(abs(recv - expected))
  call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  if (rank == 0) then
    if (global_err <= 1.0d-12) then
      write(*,'(a,i0,a,i0,a,es12.4)') &
        'cuda_aware_mpi_probe=pass ranks=', nranks, ' n=', n, ' err=', global_err
    else
      write(*,'(a,i0,a,i0,a,es12.4)') &
        'cuda_aware_mpi_probe=fail ranks=', nranks, ' n=', n, ' err=', global_err
    end if
  end if

  deallocate(recv, send)
  call MPI_Finalize(ierr)
  if (global_err > 1.0d-12) stop 2
end program cuda_aware_mpi_probe
EOF

"${MPIFORT}" ${FCFLAGS} -o "${exe}" "${src}"
# shellcheck disable=SC2086
"${MPIRUN}" ${MPI_ARGS} -np "${RANKS}" "${exe}" "${N}"
