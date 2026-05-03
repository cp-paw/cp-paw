#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
GPU_EXPLORATION_ROOT=${GPU_EXPLORATION_ROOT:-"${HERE}/runs/${TEST}-gpu-exploration-$(date +%Y%m%d-%H%M%S)"}
NSTEPS=${NSTEPS:-1}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
EMPTY_BANDS=${EMPTY_BANDS:-128}
GPU_CASES=${GPU_CASES:-"gpu gpu_resident gpu_resident_nosync gpu_force_all gpu_3dfft gpu_managed gpu_unified gpu_off gpu_resident_off"}
CPU_CASES=${CPU_CASES:-"cpu nvpl"}

mkdir -p "${GPU_EXPLORATION_ROOT}"
echo "${GPU_EXPLORATION_ROOT}" > "${HERE}/runs/latest_gpu_exploration"

COMBINED="${GPU_EXPLORATION_ROOT}/combined_benchmark.tsv"
LOG="${GPU_EXPLORATION_ROOT}/gpu_exploration.log"
: > "${LOG}"

log() {
  printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG}"
}

append_suite() {
  local suite=$1
  local tsv=$2
  [[ -f "${tsv}" ]] || return 0
  if [[ ! -s "${COMBINED}" ]]; then
    awk 'NR == 1 { print "suite\t" $0; next } { print suite "\t" $0 }' \
      suite="${suite}" "${tsv}" > "${COMBINED}"
  else
    awk 'NR > 1 { print suite "\t" $0 }' suite="${suite}" "${tsv}" >> "${COMBINED}"
  fi
}

run_suite() {
  local suite=$1
  local ranks=$2
  local cases=$3
  local root="${GPU_EXPLORATION_ROOT}/${suite}"
  local suite_log="${GPU_EXPLORATION_ROOT}/${suite}.log"

  log "START suite=${suite} ranks=${ranks} cases=${cases}"
  if env TEST="${TEST}" EMPTY_BANDS="${EMPTY_BANDS}" NSTEPS="${NSTEPS}" \
      RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${cases}" TIMEOUT="${TIMEOUT}" \
      RUN_ROOT="${root}" "${HERE}/run_benchmark.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${root}.status"
  fi
  append_suite "${suite}" "${root}/benchmark.tsv"
}

if [[ -x "${HERE}/../../../src/Tools/Scripts/paw_gpu_capabilities.sh" ]]; then
  "${HERE}/../../../src/Tools/Scripts/paw_gpu_capabilities.sh" \
    > "${GPU_EXPLORATION_ROOT}/gpu_capabilities.txt" 2>&1 || true
fi

run_suite "one_rank_gpu" 1 "${GPU_CASES}"
run_suite "one_rank_cpu" 1 "${CPU_CASES}"
run_suite "eight_rank_cpu" 8 "${CPU_CASES}"

log "ALL DONE root=${GPU_EXPLORATION_ROOT}"
[[ -f "${COMBINED}" ]] && log "combined=${COMBINED}"
