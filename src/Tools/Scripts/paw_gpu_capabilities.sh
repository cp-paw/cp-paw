#!/usr/bin/env bash
set -euo pipefail

nvhpc_platform() {
  case "$(uname -s)_$(uname -m)" in
    Linux_aarch64|Linux_arm64) echo "Linux_aarch64" ;;
    Linux_x86_64) echo "Linux_x86_64" ;;
    *) echo "" ;;
  esac
}

find_nvhpc_root() {
  local platform compiler root
  if [[ -n "${NVHPC_ROOT:-}" && -d "${NVHPC_ROOT}" ]]; then
    echo "${NVHPC_ROOT}"
    return 0
  fi
  compiler=$(command -v nvfortran 2>/dev/null || true)
  if [[ -n "${compiler}" ]]; then
    dirname "$(dirname "${compiler}")"
    return 0
  fi
  platform=$(nvhpc_platform)
  for root in /opt/nvidia/hpc_sdk/${platform}/*; do
    [[ -x "${root}/compilers/bin/nvfortran" ]] || continue
    echo "${root}"
  done | sort -V | tail -1
}

find_lib() {
  local lib=$1
  local root=${2:-}
  local dir
  for dir in \
      "${root}"/math_libs/*/targets/*/lib \
      "${root}"/REDIST/math_libs/*/targets/*/lib \
      "${root}"/math_libs/*/lib64 \
      "${root}"/comm_libs/*/nccl/lib \
      "${root}"/comm_libs/*/nvshmem/lib \
      /usr/local/cuda/lib64 \
      /usr/local/cuda-*/lib64; do
    [[ -f "${dir}/${lib}" ]] || continue
    echo "${dir}/${lib}"
    return 0
  done
}

yesno_path() {
  local key=$1
  local path=$2
  if [[ -n "${path}" ]]; then
    printf '%s=yes path=%s\n' "${key}" "${path}"
  else
    printf '%s=no\n' "${key}"
  fi
}

cuda_aware_mpi() {
  local mpirun=${MPIRUN:-}
  local ompi_info
  if [[ -z "${mpirun}" ]]; then
    mpirun=$(command -v mpirun 2>/dev/null || true)
  fi
  if [[ -z "${mpirun}" ]]; then
    echo "cuda_aware_mpi=unknown reason=no_mpirun"
    return 0
  fi
  ompi_info=$(dirname "${mpirun}")/ompi_info
  if [[ ! -x "${ompi_info}" ]]; then
    ompi_info=$(command -v ompi_info 2>/dev/null || true)
  fi
  if [[ -z "${ompi_info}" ]]; then
    echo "cuda_aware_mpi=unknown reason=no_ompi_info"
    return 0
  fi
  if "${ompi_info}" --parsable --all 2>/dev/null \
      | grep -Eiq 'cuda_support:value:true|built_with_cuda_support:value:true'; then
    echo "cuda_aware_mpi=yes mpirun=${mpirun}"
  else
    echo "cuda_aware_mpi=unknown mpirun=${mpirun}"
  fi
}

root=$(find_nvhpc_root || true)
echo "date=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "hostname=$(hostname)"
echo "os=$(uname -s)"
echo "arch=$(uname -m)"
echo "nvhpc_root=${root}"
echo "nvfortran=$(command -v nvfortran 2>/dev/null || true)"
if command -v nvfortran >/dev/null 2>&1; then
  nvfortran --version 2>&1 | head -3 | sed 's/^/nvfortran_version=/'
fi
echo "nvcc=$(command -v nvcc 2>/dev/null || true)"
if command -v nvcc >/dev/null 2>&1; then
  nvcc --version 2>&1 | tail -1 | sed 's/^/nvcc_version=/'
fi
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=index,name,compute_cap,memory.total \
    --format=csv,noheader 2>/dev/null | sed 's/^/gpu=/'
else
  echo "gpu=none"
fi

yesno_path "cublas" "$(find_lib libcublas.so "${root}" || true)"
yesno_path "cublaslt" "$(find_lib libcublasLt.so "${root}" || true)"
yesno_path "cufft" "$(find_lib libcufft.so "${root}" || true)"
yesno_path "cufftw" "$(find_lib libcufftw.so "${root}" || true)"
yesno_path "cusolver" "$(find_lib libcusolver.so "${root}" || true)"
yesno_path "cusparse" "$(find_lib libcusparse.so "${root}" || true)"
yesno_path "cutensor" "$(find_lib libcutensor.so "${root}" || true)"
yesno_path "cudss" "$(find_lib libcudss.so "${root}" || true)"
yesno_path "nccl" "$(find_lib libnccl.so "${root}" || true)"
yesno_path "nvshmem" "$(find_lib libnvshmem_host.so "${root}" || true)"
cuda_aware_mpi

echo "recommended_gpu_cases=gpu gpu_force_all gpu_3dfft gpu_off"
echo "recommended_resource_cases=cpu nvpl gpu gpu_off"
