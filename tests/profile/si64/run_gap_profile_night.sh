#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

export TEST=${TEST:-si64_bands}
export EMPTY_BANDS=${EMPTY_BANDS:-2048}
export NSTEPS=${NSTEPS:-3}
export REPEATS=${REPEATS:-1}
export TIMEOUT=${TIMEOUT:-43200}
export GPU_RANKS=${GPU_RANKS:-1}
export CPU_RANKS=${CPU_RANKS:-8}
export RUN_GPU_ALL=${RUN_GPU_ALL:-no}
export GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_addpro_host gpu_resident_pro_host gpu_off"}
export AUTO_BUILD_TARGETS=${AUTO_BUILD_TARGETS:-no}
export AUTO_BUILD_JOBS=${AUTO_BUILD_JOBS:-16}

RUN_CPU_REFERENCES=${RUN_CPU_REFERENCES:-no}
if [[ ${RUN_CPU_REFERENCES} == yes || ${RUN_CPU_REFERENCES} == true || ${RUN_CPU_REFERENCES} == 1 ]]; then
  export CPU_CASES=${CPU_CASES:-"cpu nvhpc_cpu"}
else
  export CPU_CASES=${CPU_CASES:-}
fi

export NVHPC_STANDARD_ROOT=${NVHPC_STANDARD_ROOT:-"${HERE}/runs/${TEST}-gap-profile-night-$(date +%Y%m%d-%H%M%S)"}

mkdir -p "${HERE}/runs"
echo "${NVHPC_STANDARD_ROOT}" > "${HERE}/runs/latest_gap_profile_night"

"${HERE}/run_nvhpc_standard.sh"
