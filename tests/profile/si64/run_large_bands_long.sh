#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

export TEST=${TEST:-si64_bands}
export FOLLOWUP_ROOT=${LARGE_BANDS_ROOT:-"${HERE}/runs/${TEST}-large-bands-long-$(date +%Y%m%d-%H%M%S)"}
export NSTEPS_LIST=${NSTEPS_LIST:-"3 10"}
export EMPTY_BANDS_LIST=${EMPTY_BANDS_LIST:-"512 1024"}
export REPEATS=${REPEATS:-1}
export GPU_RANKS=${GPU_RANKS:-1}
export CPU_RANKS=${CPU_RANKS:-8}
export GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_addpro_host gpu_resident_pro_host gpu_resident_nosync gpu_off"}
export ONE_RANK_CPU_CASES=${ONE_RANK_CPU_CASES:-"cpu nvhpc_cpu"}
export CPU_CASES=${CPU_CASES:-"cpu nvhpc_cpu"}
export TIMEOUT=${TIMEOUT:-14400}

"${HERE}/run_followup.sh"
