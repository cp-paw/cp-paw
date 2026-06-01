#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
. "${HERE}/case_recommendations.sh"

TEST=${TEST:-si64_bands}
EMPTY_BANDS=${EMPTY_BANDS:-2048}
NSTEPS=${NSTEPS:-1}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
GPU_RANKS=${GPU_RANKS:-1}
CPU_RANKS=${CPU_RANKS:-8}
RUN_GPU_ALL=${RUN_GPU_ALL:-no}
ADD_GPU_FALLBACK=${ADD_GPU_FALLBACK:-no}
AUTO_BUILD_TARGETS=${AUTO_BUILD_TARGETS:-no}
AUTO_BUILD_JOBS=${AUTO_BUILD_JOBS:-16}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-16}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-16}

DEFAULT_GPU_CASES="gpu_resident_stack gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_force_dedpro gpu_resident_stack_serial3dfft_accmap_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache"
if [[ -z "${GPU_CASES+x}" ]]; then
  GPU_CASES=auto
fi
if [[ -z "${CPU_CASES+x}" ]]; then
  CPU_CASES=auto
fi

GPU_CASES=$(cppaw_resolve_recommended_cases \
  "${GPU_CASES}" recommended_large_band_gpu_cases "${DEFAULT_GPU_CASES}")
CPU_CASES=$(cppaw_resolve_recommended_cases \
  "${CPU_CASES}" recommended_cpu_cases "cpu nvhpc_cpu")

export TEST EMPTY_BANDS NSTEPS REPEATS TIMEOUT
export GPU_RANKS CPU_RANKS GPU_CASES CPU_CASES
export RUN_GPU_ALL ADD_GPU_FALLBACK AUTO_BUILD_TARGETS AUTO_BUILD_JOBS
export PROFILE_ROW_TOP PRESENT_ROW_TOP
export NVHPC_STANDARD_ROOT=${GPU_RESOURCE_COMPARISON_ROOT:-"${HERE}/runs/${TEST}-gpu-resource-comparison-$(date +%Y%m%d-%H%M%S)"}

mkdir -p "${NVHPC_STANDARD_ROOT}"
echo "${NVHPC_STANDARD_ROOT}" > "${HERE}/runs/latest_gpu_resource_comparison"

exec "${HERE}/run_nvhpc_standard.sh"
