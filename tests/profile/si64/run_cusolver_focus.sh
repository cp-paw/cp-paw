#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
CUSOLVER_FOCUS_ROOT=${CUSOLVER_FOCUS_ROOT:-"${HERE}/runs/${TEST}-cusolver-focus-$(date +%Y%m%d-%H%M%S)"}
NSTEPS=${NSTEPS:-1}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
EMPTY_BANDS_LIST=${EMPTY_BANDS_LIST:-"128 256 512"}
GPU_RANKS=${GPU_RANKS:-1}
CPU_RANKS=${CPU_RANKS:-8}
CUSOLVER_CASES=${CUSOLVER_CASES:-"cusolver cusolver_generalized cusolver_generalized_conservative cusolver_off"}
CPU_CASES=${CPU_CASES:-"cpu nvhpc_cpu"}

mkdir -p "${CUSOLVER_FOCUS_ROOT}"
echo "${CUSOLVER_FOCUS_ROOT}" > "${HERE}/runs/latest_cusolver_focus"

COMBINED="${CUSOLVER_FOCUS_ROOT}/combined_benchmark.tsv"
LOG="${CUSOLVER_FOCUS_ROOT}/cusolver_focus.log"
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
  local empty_bands=$2
  local ranks=$3
  local cases=$4
  local root="${CUSOLVER_FOCUS_ROOT}/${suite}"
  local suite_log="${CUSOLVER_FOCUS_ROOT}/${suite}.log"

  log "START suite=${suite} empty_bands=${empty_bands} ranks=${ranks} cases=${cases}"
  if env TEST="${TEST}" EMPTY_BANDS="${empty_bands}" NSTEPS="${NSTEPS}" \
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

for empty_bands in ${EMPTY_BANDS_LIST}; do
  run_suite "empty${empty_bands}_${GPU_RANKS}rank_cusolver" \
    "${empty_bands}" "${GPU_RANKS}" "${CUSOLVER_CASES}"
  run_suite "empty${empty_bands}_1rank_cpu" \
    "${empty_bands}" 1 "${CPU_CASES}"
  run_suite "empty${empty_bands}_${CPU_RANKS}rank_cpu_ref" \
    "${empty_bands}" "${CPU_RANKS}" "${CPU_CASES}"
done

log "ALL DONE root=${CUSOLVER_FOCUS_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${CUSOLVER_FOCUS_ROOT}/combined_benchmark.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${CUSOLVER_FOCUS_ROOT}/combined_compare.md" || true
  log "combined=${COMBINED}"
fi
