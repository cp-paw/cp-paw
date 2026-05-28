#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
NVHPC_STANDARD_ROOT=${NVHPC_STANDARD_ROOT:-"${HERE}/runs/${TEST}-nvhpc-standard-$(date +%Y%m%d-%H%M%S)"}
NSTEPS=${NSTEPS:-3}
EMPTY_BANDS=${EMPTY_BANDS:-1024}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
GPU_RANKS=${GPU_RANKS:-1}
CPU_RANKS=${CPU_RANKS:-8}
GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_invbatch_off gpu_resident_no_cusolver gpu_off"}
CPU_CASES=${CPU_CASES:-"cpu nvhpc_cpu"}

mkdir -p "${NVHPC_STANDARD_ROOT}"
echo "${NVHPC_STANDARD_ROOT}" > "${HERE}/runs/latest_nvhpc_standard"

COMBINED="${NVHPC_STANDARD_ROOT}/combined_benchmark.tsv"
LOG="${NVHPC_STANDARD_ROOT}/nvhpc_standard.log"
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
  local root="${NVHPC_STANDARD_ROOT}/${suite}"
  local suite_log="${NVHPC_STANDARD_ROOT}/${suite}.log"

  log "START suite=${suite} empty_bands=${EMPTY_BANDS} nsteps=${NSTEPS} ranks=${ranks} cases=${cases}"
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

run_suite "gpu_${GPU_RANKS}rank" "${GPU_RANKS}" "${GPU_CASES}"
run_suite "cpu_1rank" 1 "${CPU_CASES}"
run_suite "cpu_${CPU_RANKS}rank_ref" "${CPU_RANKS}" "${CPU_CASES}"

log "ALL DONE root=${NVHPC_STANDARD_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${NVHPC_STANDARD_ROOT}/combined_benchmark.md" || true
  log "combined=${COMBINED}"
fi
