#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
TEST=${TEST:-si64}
CNTL_FILE=${CNTL_FILE:-"${HERE}/${TEST}.cntl"}
STRC_FILE=${STRC_FILE:-"${HERE}/${TEST}.strc"}
if [[ ! -f "${STRC_FILE}" && -f "${HERE}/si64.strc" ]]; then
  STRC_FILE="${HERE}/si64.strc"
fi
NSTEPS=${NSTEPS:-1}
RANKS=${RANKS:-1}
CASE=${CASE:-gpu_resident}
TIMEOUT=${TIMEOUT:-600}
RUN_ROOT=${RUN_ROOT:-"${HERE}/runs/${TEST}-${CASE}-nsys-nstep${NSTEPS}-${RANKS}ranks-$(date +%Y%m%d-%H%M%S)"}
MPI_ARGS=${MPI_ARGS:---mca coll ^hcoll}
NSYS_TRACE=${NSYS_TRACE:-cuda,nvtx,osrt}
NSYS_STATS=${NSYS_STATS:-true}
NSYS_MPI_MODE=${NSYS_MPI_MODE:-outer}
TIMEOUT_PREFIX=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_PREFIX="timeout ${TIMEOUT}s"
fi

nvhpc_platform() {
  case "$(uname -s)_$(uname -m)" in
    Linux_aarch64|Linux_arm64) echo "Linux_aarch64" ;;
    Linux_x86_64) echo "Linux_x86_64" ;;
    *) echo "" ;;
  esac
}

default_mpirun() {
  local platform root candidate
  platform=$(nvhpc_platform)
  for root in "${NVHPC_ROOT:-}" /opt/nvidia/hpc_sdk/${platform}/*; do
    [[ -n "${root}" && -d "${root}" ]] || continue
    for candidate in \
        "${root}/comm_libs/hpcx/bin/mpirun" \
        "${root}"/comm_libs/*/hpcx/*/ompi/bin/mpirun \
        "${root}"/comm_libs/*/hpcx/bin/mpirun; do
      if [[ -x "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    done
  done
  command -v mpirun 2>/dev/null || echo mpirun
}

MPIRUN=${MPIRUN:-$(default_mpirun)}

serial_exe() {
  case "$1" in
    gpu_all*) echo "${ROOT}/bin/nvhpc_gpu_all_profile/paw_nvhpc_gpu_all_profile.x" ;;
    gpu_resident*) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" ;;
    gpu_managed*) echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile/paw_nvhpc_gpu_acc_managed_profile.x" ;;
    gpu_unified*) echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile/paw_nvhpc_gpu_acc_unified_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile/paw_nvhpc_cublas_acc_profile.x" ;;
    cusolver*) echo "${ROOT}/bin/nvhpc_cusolver_acc_profile/paw_nvhpc_cusolver_acc_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile/paw_nvhpc_cufftw_profile.x" ;;
    cufft*) echo "${ROOT}/bin/nvhpc_cufft_profile/paw_nvhpc_cufft_profile.x" ;;
    nvlamath) echo "${ROOT}/bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" ;;
    nvhpc_cpu|nvpl) echo "${ROOT}/bin/nvhpc_profile/paw_nvhpc_profile.x" ;;
    *) echo "unknown Nsight case ${1}" >&2; return 1 ;;
  esac
}

parallel_exe() {
  case "$1" in
    gpu_all*) echo "${ROOT}/bin/nvhpc_gpu_all_profile_parallel/ppaw_nvhpc_gpu_all_profile.x" ;;
    gpu_resident*) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile_parallel/ppaw_nvhpc_gpu_acc_residency_profile.x" ;;
    gpu_managed*) echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile_parallel/ppaw_nvhpc_gpu_acc_managed_profile.x" ;;
    gpu_unified*) echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile_parallel/ppaw_nvhpc_gpu_acc_unified_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_gpu_acc_profile_parallel/ppaw_nvhpc_gpu_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" ;;
    cusolver*) echo "${ROOT}/bin/nvhpc_cusolver_acc_profile_parallel/ppaw_nvhpc_cusolver_acc_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" ;;
    cufft*) echo "${ROOT}/bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" ;;
    nvlamath) echo "${ROOT}/bin/nvhpc_nvlamath_profile_parallel/ppaw_nvhpc_nvlamath_profile.x" ;;
    nvhpc_cpu|nvpl) echo "${ROOT}/bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" ;;
    *) echo "unknown Nsight case ${1}" >&2; return 1 ;;
  esac
}

case_exe() {
  if [[ "${RANKS}" -gt 1 ]]; then
    parallel_exe "$1"
  else
    serial_exe "$1"
  fi
}

cufft_min_default() {
  echo "${CPPAW_CUFFT_CONSERVATIVE_MIN_ELEMENTS:-1000000}"
}

cufft_env() {
  local min_elements
  min_elements=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-$(cufft_min_default)}
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${min_elements}"
}

cufft_force_env() {
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0}"
}

cufft3d_env() {
  local min_elements min_3d_elements
  min_elements=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-$(cufft_min_default)}
  min_3d_elements=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-$(cufft_min_default)}
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${min_elements} CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${min_3d_elements}"
}

cublas_env() {
  echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_ACC_MINFLOP:-1e7}"
}

cublas_conservative_env() {
  echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}"
}

cublas_projection_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_PROJECTION_MINFLOP=${CPPAW_CUBLAS_PROJECTION_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_overlap_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_OVERLAP_MINFLOP=${CPPAW_CUBLAS_OVERLAP_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_addproduct_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP=${CPPAW_CUBLAS_ADDPRODUCT_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_matmul_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_MATMUL_MINFLOP=${CPPAW_CUBLAS_MATMUL_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cusolver_env() {
  local min_n=$1
  local env_line="CPPAW_CUSOLVER_ACC_MIN_N=${min_n}"
  if [[ -n "${CPPAW_CUSOLVER_ACC_STANDARD_MIN_N:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_ACC_STANDARD_MIN_N}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK=${CPPAW_CUSOLVER_ACC_CHECK}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK_TOL:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK_TOL=${CPPAW_CUSOLVER_ACC_CHECK_TOL}"
  fi
  echo "${env_line}"
}

cusolver_standard_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_FORCE_MIN_N:-1}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_OFF_MIN_N:-1000000000}"
  echo "${env_line}"
}

cusolver_generalized_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_OFF_MIN_N:-1000000000}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_FORCE_MIN_N:-1}"
  echo "${env_line}"
}

cusolver_generalized_conservative_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_OFF_MIN_N:-1000000000}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_CONSERVATIVE_MIN_N:-256}"
  echo "${env_line}"
}

case_env() {
  case "$1" in
    cublas) cublas_env ;;
    cublas_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    cublas_invbatch_off) echo "$(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0" ;;
    cublas_conservative) cublas_conservative_env ;;
    cublas_projection_conservative) cublas_projection_conservative_env ;;
    cublas_overlap_conservative) cublas_overlap_conservative_env ;;
    cublas_addproduct_conservative) cublas_addproduct_conservative_env ;;
    cublas_matmul_conservative) cublas_matmul_conservative_env ;;
    cublas_off) echo "CPPAW_CUBLAS_ACC=0" ;;
    cusolver) cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}" ;;
    cusolver_standard) cusolver_standard_env ;;
    cusolver_generalized) cusolver_generalized_env ;;
    cusolver_generalized_conservative) cusolver_generalized_conservative_env ;;
    cusolver_conservative) cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}" ;;
    cusolver_off) echo "CPPAW_CUSOLVER_ACC=0" ;;
    cufft) cufft_env ;;
    cufft_force_all) cufft_force_env ;;
    cufft_off) echo "CPPAW_CUFFT_ACC=0" ;;
    gpu) cublas_env ;;
    gpu_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_invbatch_off) echo "$(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0" ;;
    gpu_projection_conservative) cublas_projection_conservative_env ;;
    gpu_overlap_conservative) cublas_overlap_conservative_env ;;
    gpu_addproduct_conservative) cublas_addproduct_conservative_env ;;
    gpu_matmul_conservative) cublas_matmul_conservative_env ;;
    gpu_force_all) echo "$(cufft_force_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_conservative) echo "$(cublas_conservative_env) $(cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}")" ;;
    gpu_resident) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_nosync) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_resident_orthox) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_orthox_off) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_orthox_nosync) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_resident_orthoconst) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_CONST_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_opsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_denmat_energy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_denmat_energy_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_batch) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_batch) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_hpsi_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_psim_propagate|gpu_resident_psim) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 $(cublas_env)" ;;
    gpu_hpsi_psim_propagate|gpu_resident_hpsi_psim) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 $(cublas_env)" ;;
    gpu_resident_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cublas_env)" ;;
    gpu_resident_stack_psi0_ortho_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_ORTHO_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_psi0_prinfo_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_PRINFO_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_setup_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_SETUP_PSI_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_psim_phase) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_hpsi_prop_psim_phase) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_hpsi_prop_psim_switch) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1 CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_cufft) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cufft_env) $(cublas_env)" ;;
    gpu_resident_stack_cufft_force) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cufft_force_env) $(cublas_env)" ;;
    gpu_resident_addpro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE=0 $(cublas_env)" ;;
    gpu_resident_addpro_hpsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_HPSI=0 $(cublas_env)" ;;
    gpu_resident_addpro_opsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_OPSI=0 $(cublas_env)" ;;
    gpu_resident_opsi_addpro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_OPSI=0 $(cublas_env)" ;;
    gpu_resident_pro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PRO_EXPANSION=0 $(cublas_env)" ;;
    gpu_resident_forcepsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_FORCE_PSI_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_invbatch_off) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_CUBLAS_ACC_INVERSION_BATCH=0 $(cublas_env)" ;;
    gpu_resident_no_cusolver) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUSOLVER_ACC=0" ;;
    gpu_resident_1coverlap) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_1COVERLAP=1 $(cublas_env)" ;;
    gpu_resident_1coverlap_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_1COVERLAP=0 $(cublas_env)" ;;
    gpu_resident_gram_cholesky) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GRAM_CHOLESKY=1 $(cublas_env)" ;;
    gpu_resident_gram_legacy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GRAM_CHOLESKY=0 $(cublas_env)" ;;
    gpu_resident_projection_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_projection_conservative_env)" ;;
    gpu_resident_overlap_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_overlap_conservative_env)" ;;
    gpu_resident_addproduct_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_addproduct_conservative_env)" ;;
    gpu_resident_matmul_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_matmul_conservative_env)" ;;
    gpu_resident_force_all) echo "CPPAW_GPU_RESIDENCY=1 $(cufft_force_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_resident_off) echo "CPPAW_GPU_RESIDENCY=0 CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_all) echo "$(cufft_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_nosync) echo "$(cufft_env) $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_invbatch_off) echo "$(cufft_env) $(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_no_cufft) echo "CPPAW_CUFFT_ACC=0 $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_no_cublas) echo "$(cufft_force_env) CPPAW_CUBLAS_ACC=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_no_cusolver) echo "$(cufft_force_env) $(cublas_env) CPPAW_CUSOLVER_ACC=0" ;;
    gpu_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_managed|gpu_unified) cublas_env ;;
    nvlamath|cufftw|nvhpc_cpu|nvpl) echo "" ;;
    *) echo "unknown Nsight case ${1}" >&2; return 1 ;;
  esac
}

EXE=${EXE:-$(case_exe "${CASE}")}
CASE_ENV=${CASE_ENV:-$(case_env "${CASE}")}

if ! command -v nsys >/dev/null 2>&1; then
  echo "nsys was not found in PATH." >&2
  exit 1
fi
if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${EXE}" >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}"
if [[ ! -f "${CNTL_FILE}" ]]; then
  echo "Control file not found: ${CNTL_FILE}" >&2
  exit 1
fi
if [[ ! -f "${STRC_FILE}" ]]; then
  echo "Structure file not found: ${STRC_FILE}" >&2
  exit 1
fi
cp "${CNTL_FILE}" "${RUN_ROOT}/${TEST}.cntl"
cp "${STRC_FILE}" "${RUN_ROOT}/${TEST}.strc"
cp "${HERE}/profile_summary.py" "${RUN_ROOT}/"
cp "${HERE}/nsys_sql_summary.py" "${RUN_ROOT}/"
cp "${ROOT}/tests/fulltests/si2/stp.cntl" "${RUN_ROOT}/"
perl -0pi -e "s/NSTEP\\s*=\\s*\\d+/NSTEP=${NSTEPS}/" "${RUN_ROOT}/${TEST}.cntl"
if [[ -n "${EMPTY_BANDS:-}" ]]; then
  perl -0pi -e "s/EMPTY\\s*=\\s*\\d+/EMPTY=${EMPTY_BANDS}/" "${RUN_ROOT}/${TEST}.strc"
fi

cd "${RUN_ROOT}"
echo "${RUN_ROOT}" > "${HERE}/runs/latest_nsys"
{
  echo "CASE=${CASE}"
  echo "CASE_ENV=${CASE_ENV}"
  echo "EXE=${EXE}"
  echo "RANKS=${RANKS}"
  echo "NSYS_TRACE=${NSYS_TRACE}"
} > nsys_case.env

if [[ "${RANKS}" -gt 1 && "${NSYS_MPI_MODE}" == "per_rank" ]]; then
  cat > nsys_rank_wrapper.sh <<EOF
#!/usr/bin/env bash
set -euo pipefail
rank=\${OMPI_COMM_WORLD_RANK:-\${PMI_RANK:-\${SLURM_PROCID:-0}}}
exec env CPPAW_ACCEL_PROFILE_FILE=nsys_profile ${CASE_ENV} \\
  nsys profile --force-overwrite=true --stats="${NSYS_STATS}" --trace="${NSYS_TRACE}" \\
  -o "nsys_rank\${rank}" "${EXE}" "${TEST}.cntl"
EOF
  chmod +x nsys_rank_wrapper.sh
  # shellcheck disable=SC2086
  CMD="${MPIRUN} ${MPI_ARGS} -np ${RANKS} ./nsys_rank_wrapper.sh"
elif [[ "${RANKS}" -gt 1 ]]; then
  # CP-PAW exits MPI runs via MPI_ABORT(0), so keep nsys outside the MPI job
  # by default to let it flush the combined report.
  CMD="env CPPAW_ACCEL_PROFILE_FILE=nsys_profile ${CASE_ENV} nsys profile --force-overwrite=true --stats=${NSYS_STATS} --trace=${NSYS_TRACE} -o nsys_mpi ${MPIRUN} ${MPI_ARGS} -np ${RANKS} ${EXE} ${TEST}.cntl"
else
  CMD="env CPPAW_ACCEL_PROFILE_FILE=nsys_profile ${CASE_ENV} nsys profile --force-overwrite=true --stats=${NSYS_STATS} --trace=${NSYS_TRACE} -o nsys ${EXE} ${TEST}.cntl"
fi

echo "running Nsight Systems case ${CASE}: ${CMD}"
# shellcheck disable=SC2086
/usr/bin/time -p ${TIMEOUT_PREFIX} ${CMD} > out.log 2> err.log

python3 profile_summary.py nsys_profile*.csv > summary.txt
if compgen -G "*.sqlite" >/dev/null; then
  python3 nsys_sql_summary.py ./*.sqlite > nsys_sql_summary.txt || true
fi
grep -E "^real |^user |^sys " err.log > time.txt || true
echo "Nsight data: ${RUN_ROOT}"
