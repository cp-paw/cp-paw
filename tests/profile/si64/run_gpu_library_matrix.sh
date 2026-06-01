#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

export TEST=${TEST:-si64_bands}
export EMPTY_BANDS=${EMPTY_BANDS:-1024}
export NSTEPS=${NSTEPS:-1}
export REPEATS=${REPEATS:-1}
export TIMEOUT=${TIMEOUT:-7200}
export GPU_RANKS=${GPU_RANKS:-1}
export CPU_RANKS=${CPU_RANKS:-8}
export RUN_GPU_ALL=${RUN_GPU_ALL:-no}
export GPU_CASES=${GPU_CASES:-"gpu_off gpu_all gpu_all_off gpu_resident gpu_resident_orthox_off gpu_resident_no_cusolver cublas cusolver cufft cufftw nvlamath nvblas gpu_no_cufft gpu_no_cublas gpu_no_cusolver gpu_managed gpu_unified"}
export CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}
export AUTO_BUILD_TARGETS=${AUTO_BUILD_TARGETS:-no}
export AUTO_BUILD_JOBS=${AUTO_BUILD_JOBS:-16}
export NVHPC_STANDARD_ROOT=${NVHPC_STANDARD_ROOT:-${GPU_LIBRARY_MATRIX_ROOT:-"${HERE}/runs/${TEST}-gpu-library-matrix-$(date +%Y%m%d-%H%M%S)"}}

mkdir -p "${HERE}/runs"
echo "${NVHPC_STANDARD_ROOT}" > "${HERE}/runs/latest_gpu_library_matrix"

"${HERE}/run_nvhpc_standard.sh"
