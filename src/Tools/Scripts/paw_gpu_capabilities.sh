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

find_nvhpc_mpirun() {
  local root=${1:-}
  local platform candidate version_root
  platform=$(nvhpc_platform)

  for candidate in \
      "${MPIRUN:-}" \
      "${root}"/comm_libs/hpcx/bin/mpirun \
      "${root}"/comm_libs/*/hpcx/*/ompi/bin/mpirun \
      "${root}"/comm_libs/*/hpcx/bin/mpirun \
      "${root}"/../*/comm_libs/hpcx/bin/mpirun \
      "${root}"/../*/comm_libs/*/hpcx/*/ompi/bin/mpirun \
      "${root}"/../*/comm_libs/*/hpcx/bin/mpirun \
      /opt/nvidia/hpc_sdk/${platform}/*/comm_libs/hpcx/bin/mpirun \
      /opt/nvidia/hpc_sdk/${platform}/*/comm_libs/*/hpcx/*/ompi/bin/mpirun \
      /opt/nvidia/hpc_sdk/${platform}/*/comm_libs/*/hpcx/bin/mpirun \
      "${HOME:-}"/opt/nvidia/hpc_sdk/${platform}/*/comm_libs/hpcx/bin/mpirun \
      "${HOME:-}"/opt/nvidia/hpc_sdk/${platform}/*/comm_libs/*/hpcx/*/ompi/bin/mpirun \
      "${HOME:-}"/opt/nvidia/hpc_sdk/${platform}/*/comm_libs/*/hpcx/bin/mpirun; do
    if [[ -x "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done

  version_root=$(command -v mpirun 2>/dev/null || true)
  [[ -n "${version_root}" ]] && echo "${version_root}"
}

find_host_fftw_include() {
  local root=${1:-}
  local dir
  local pc_includedir

  if [[ -n "${CPPAW_FFTW3_INCLUDE:-}" && -f "${CPPAW_FFTW3_INCLUDE}" ]]; then
    echo "${CPPAW_FFTW3_INCLUDE}"
    return 0
  fi

  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists fftw3; then
    pc_includedir=$(pkg-config --variable=includedir fftw3 2>/dev/null || true)
    if [[ -f "${pc_includedir}/fftw3.f03" ]]; then
      echo "${pc_includedir}/fftw3.f03"
      return 0
    fi
  fi

  for dir in \
      "${root}"/compilers/include/nvpl_fftw \
      "${root}"/math_libs/nvpl/include/nvpl_fftw \
      "${root}"/REDIST/math_libs/nvpl/include/nvpl_fftw \
      /usr/local/include \
      /usr/include \
      /opt/homebrew/include; do
    if [[ -f "${dir}/fftw3.f03" ]]; then
      echo "${dir}/fftw3.f03"
      return 0
    fi
  done
}

find_host_fftw() {
  local root=${1:-}
  local dir
  local include_path
  local pc_prefix

  include_path=$(find_host_fftw_include "${root}" || true)
  [[ -n "${include_path}" ]] || return 1

  if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists fftw3; then
    pc_prefix=$(pkg-config --variable=prefix fftw3 2>/dev/null || true)
    echo "pkg-config:fftw3${pc_prefix:+ prefix=${pc_prefix}} include=${include_path}"
    return 0
  fi

  for dir in \
      "${root}"/math_libs/nvpl/lib \
      "${root}"/REDIST/math_libs/nvpl/lib \
      "${root}"/compilers/lib \
      "${root}"/REDIST/compilers/lib \
      /usr/local/lib \
      /usr/lib64 \
      /usr/lib \
      /opt/homebrew/lib; do
    if [[ -f "${dir}/libnvpl_fftw.so" ]]; then
      echo "${dir}/libnvpl_fftw.so include=${include_path}"
      return 0
    fi
    if [[ -f "${dir}/libfftw3.so" || -f "${dir}/libfftw3.dylib" ]]; then
      echo "${dir}/libfftw3 include=${include_path}"
      return 0
    fi
  done
}

find_host_blas_lapack() {
  local root=${1:-}
  local dir
  local pc
  local pc_prefix

  if [[ "$(uname -s)" == Darwin ]]; then
    echo "framework:Accelerate"
    return 0
  fi

  if command -v pkg-config >/dev/null 2>&1; then
    for pc in openblas mkl; do
      if pkg-config --exists "${pc}"; then
        pc_prefix=$(pkg-config --variable=prefix "${pc}" 2>/dev/null || true)
        echo "pkg-config:${pc}${pc_prefix:+ prefix=${pc_prefix}}"
        return 0
      fi
    done
    if pkg-config --exists lapack && pkg-config --exists blas; then
      pc_prefix=$(pkg-config --variable=prefix lapack 2>/dev/null || true)
      echo "pkg-config:lapack+blas${pc_prefix:+ prefix=${pc_prefix}}"
      return 0
    fi
  fi

  for dir in \
      "${root}"/math_libs/nvpl/lib \
      "${root}"/REDIST/math_libs/nvpl/lib \
      "${root}"/compilers/lib \
      "${root}"/REDIST/compilers/lib \
      /usr/local/lib \
      /usr/lib64 \
      /usr/lib \
      /opt/homebrew/lib; do
    if [[ -f "${dir}/libnvpl_blas_lp64_seq.so" \
          && -f "${dir}/libnvpl_lapack_lp64_seq.so" ]]; then
      echo "${dir}/libnvpl_{blas,lapack}_lp64_seq.so"
      return 0
    fi
    if [[ -f "${dir}/libblas.so" && -f "${dir}/liblapack.so" ]]; then
      echo "${dir}/lib{blas,lapack}.so"
      return 0
    fi
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
  local root=${1:-}
  local mpirun=${MPIRUN:-}
  local ompi_info
  if [[ -z "${mpirun}" ]]; then
    mpirun=$(find_nvhpc_mpirun "${root}" || true)
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
mpirun_path=$(find_nvhpc_mpirun "${root}" || true)
echo "date=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "hostname=$(hostname)"
echo "os=$(uname -s)"
echo "arch=$(uname -m)"
echo "nvhpc_root=${root}"
nvfortran_path=$(command -v nvfortran 2>/dev/null || true)
if [[ -z "${nvfortran_path}" && -n "${root}" && -x "${root}/compilers/bin/nvfortran" ]]; then
  nvfortran_path="${root}/compilers/bin/nvfortran"
fi
echo "nvfortran=${nvfortran_path}"
if [[ -n "${nvfortran_path}" ]]; then
  "${nvfortran_path}" --version 2>&1 | head -3 | sed 's/^/nvfortran_version=/'
fi
echo "mpirun=${mpirun_path}"
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
host_fftw_path=$(find_host_fftw "${root}" || true)
host_blas_lapack_path=$(find_host_blas_lapack "${root}" || true)

yesno_path "host_fftw" "${host_fftw_path}"
yesno_path "host_blas_lapack" "${host_blas_lapack_path}"
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
cuda_aware_mpi "${root}"
echo "cuda_aware_mpi_probe=run src/Tools/Scripts/paw_cuda_aware_mpi_probe.sh"

recommended_cpu_cases="cpu"
if [[ -n "${root}" && -x "${root}/compilers/bin/nvfortran" ]]; then
  recommended_cpu_cases=$(append_case "${recommended_cpu_cases}" "nvhpc_cpu")
fi
recommended_gpu_cases="none"
recommended_gpu_diagnostic_cases=""
recommended_resource_cases="${recommended_cpu_cases}"

if [[ -z "${host_fftw_path}" ]]; then
  recommended_cpu_cases="none"
  recommended_resource_cases="none"
  echo "recommended_cpu_reason=no_host_fftw_runtime_found"
  echo "recommended_gpu_reason=no_host_fftw_runtime_found"
elif [[ -z "${host_blas_lapack_path}" ]]; then
  recommended_cpu_cases="none"
  recommended_resource_cases="none"
  echo "recommended_cpu_reason=no_host_blas_lapack_runtime_found"
  echo "recommended_gpu_reason=no_host_blas_lapack_runtime_found"
elif [[ "${has_gpu}" == yes && -n "${cublas_path}" ]]; then
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
