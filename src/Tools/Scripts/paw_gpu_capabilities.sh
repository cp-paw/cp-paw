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
  local platform compiler root version_root
  platform=$(nvhpc_platform)
  if [[ -n "${NVHPC_ROOT:-}" && -x "${NVHPC_ROOT}/compilers/bin/nvfortran" ]]; then
    echo "${NVHPC_ROOT}"
    return 0
  fi
  if [[ -n "${NVHPC_ROOT:-}" && -d "${NVHPC_ROOT}" ]]; then
    version_root=$(
      for root in "${NVHPC_ROOT}"/* "${NVHPC_ROOT}/${platform}"/*; do
        [[ -x "${root}/compilers/bin/nvfortran" ]] || continue
        echo "${root}"
      done | sort -V | tail -1
    )
    if [[ -n "${version_root}" ]]; then
      echo "${version_root}"
      return 0
    fi
  fi
  compiler=$(command -v nvfortran 2>/dev/null || true)
  if [[ -n "${compiler}" ]]; then
    dirname "$(dirname "$(dirname "${compiler}")")"
    return 0
  fi
  for root in \
      /opt/nvidia/hpc_sdk/${platform}/* \
      "${HOME:-}"/opt/nvidia/hpc_sdk/${platform}/*; do
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
  "${ompi_info}" --parsable --all 2>/dev/null \
      | grep -Ei 'cuda|gpu|ucx|hcoll' \
      | sed -n '1,40s/^/mpi_capability=/' || true
}

append_case() {
  local list=$1
  local item=$2
  case " ${list} " in
    *" ${item} "*) echo "${list}" ;;
    *) echo "${list:+${list} }${item}" ;;
  esac
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
has_gpu=no
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_lines=$(nvidia-smi --query-gpu=index,name,compute_cap,memory.total \
    --format=csv,noheader 2>/dev/null || true)
  if [[ -n "${gpu_lines}" ]]; then
    has_gpu=yes
    printf '%s\n' "${gpu_lines}" | sed 's/^/gpu=/'
  else
    echo "gpu=none"
  fi
else
  echo "gpu=none"
fi

cublas_path=$(find_lib libcublas.so "${root}" || true)
cublaslt_path=$(find_lib libcublasLt.so "${root}" || true)
cufft_path=$(find_lib libcufft.so "${root}" || true)
cufftw_path=$(find_lib libcufftw.so "${root}" || true)
cusolver_path=$(find_lib libcusolver.so "${root}" || true)
cusparse_path=$(find_lib libcusparse.so "${root}" || true)
cutensor_path=$(find_lib libcutensor.so "${root}" || true)
cudss_path=$(find_lib libcudss.so "${root}" || true)
nccl_path=$(find_lib libnccl.so "${root}" || true)
nvshmem_path=$(find_lib libnvshmem_host.so "${root}" || true)

yesno_path "cublas" "${cublas_path}"
yesno_path "cublaslt" "${cublaslt_path}"
yesno_path "cufft" "${cufft_path}"
yesno_path "cufftw" "${cufftw_path}"
yesno_path "cusolver" "${cusolver_path}"
yesno_path "cusparse" "${cusparse_path}"
yesno_path "cutensor" "${cutensor_path}"
yesno_path "cudss" "${cudss_path}"
yesno_path "nccl" "${nccl_path}"
yesno_path "nvshmem" "${nvshmem_path}"
cuda_aware_mpi
echo "cuda_aware_mpi_probe=run src/Tools/Scripts/paw_cuda_aware_mpi_probe.sh"

recommended_cpu_cases="cpu nvhpc_cpu"
recommended_gpu_cases="none"
recommended_gpu_diagnostic_cases=""
recommended_resource_cases="${recommended_cpu_cases}"

if [[ "${has_gpu}" == yes && -n "${cublas_path}" ]]; then
  recommended_gpu_cases="gpu_resident_stack gpu_resident_off"
  recommended_resource_cases="${recommended_cpu_cases} gpu_resident_stack"
  recommended_gpu_diagnostic_cases="gpu_resident_stack_force_dedpro gpu_resident_nosync"
  if [[ -n "${cusolver_path}" ]]; then
    recommended_gpu_diagnostic_cases=$(append_case \
      "${recommended_gpu_diagnostic_cases}" "gpu_resident_no_cusolver")
  fi
  if [[ -n "${cufft_path}" ]]; then
    recommended_gpu_cases=$(append_case \
      "${recommended_gpu_cases}" "gpu_resident_stack_cufft")
    recommended_gpu_cases=$(append_case \
      "${recommended_gpu_cases}" "gpu_resident_stack_serial3dfft")
    recommended_gpu_diagnostic_cases=$(append_case \
      "${recommended_gpu_diagnostic_cases}" \
      "gpu_resident_stack_serial3dfft_force_dedpro")
    recommended_gpu_diagnostic_cases=$(append_case \
      "${recommended_gpu_diagnostic_cases}" \
      "gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache")
  fi
  recommended_gpu_diagnostic_cases=$(append_case \
    "${recommended_gpu_diagnostic_cases}" "gpu_force_all")
elif [[ "${has_gpu}" != yes ]]; then
  echo "recommended_gpu_reason=no_cuda_device_visible_to_nvidia_smi"
elif [[ -z "${cublas_path}" ]]; then
  echo "recommended_gpu_reason=no_cublas_runtime_found"
fi

echo "recommended_cpu_cases=${recommended_cpu_cases}"
echo "recommended_gpu_cases=${recommended_gpu_cases}"
echo "recommended_gpu_diagnostic_cases=${recommended_gpu_diagnostic_cases:-none}"
echo "recommended_resource_cases=${recommended_resource_cases}"
echo "recommended_standard_command=cd tests/profile/si64 && TEST=si64_bands EMPTY_BANDS=1024 NSTEPS=1 ./run_nvhpc_standard.sh"
