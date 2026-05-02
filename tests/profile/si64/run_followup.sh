#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
FOLLOWUP_ROOT=${FOLLOWUP_ROOT:-"${HERE}/runs/${TEST}-followup-$(date +%Y%m%d-%H%M%S)"}
NSTEPS_LIST=${NSTEPS_LIST:-"1 3 10"}
GPU_RANKS=${GPU_RANKS:-1}
CPU_RANKS=${CPU_RANKS:-8}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
EMPTY_BANDS=${EMPTY_BANDS:-128}
GPU_CASES=${GPU_CASES:-"gpu gpu_nosync gpu_off"}
CPU_CASES=${CPU_CASES:-"cpu nvpl"}
ONE_RANK_CPU_CASES=${ONE_RANK_CPU_CASES:-"cpu nvpl"}

mkdir -p "${FOLLOWUP_ROOT}"
echo "${FOLLOWUP_ROOT}" > "${HERE}/runs/latest_followup"

COMBINED="${FOLLOWUP_ROOT}/combined_benchmark.tsv"
SUMMARY_LOG="${FOLLOWUP_ROOT}/followup.log"
: > "${SUMMARY_LOG}"

log() {
  printf '%s %s\n' "$(iso_now)" "$*" | tee -a "${SUMMARY_LOG}"
}

iso_now() {
  date -u +%Y-%m-%dT%H:%M:%SZ
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
  local nsteps=$2
  local ranks=$3
  local cases=$4
  shift 4
  local suite_root="${FOLLOWUP_ROOT}/${suite}"
  local suite_log="${FOLLOWUP_ROOT}/${suite}.log"

  log "START suite=${suite} nsteps=${nsteps} ranks=${ranks} cases=${cases} env=[$*]"
  if env TEST="${TEST}" EMPTY_BANDS="${EMPTY_BANDS}" NSTEPS="${nsteps}" \
      RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${cases}" TIMEOUT="${TIMEOUT}" \
      RUN_ROOT="${suite_root}" "$@" "${HERE}/run_benchmark.sh" \
      > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${suite_root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${suite_root}.status"
  fi
  append_suite "${suite}" "${suite_root}/benchmark.tsv"
}

for nsteps in ${NSTEPS_LIST}; do
  run_suite "nstep${nsteps}_${GPU_RANKS}rank_gpu" "${nsteps}" "${GPU_RANKS}" \
    "${GPU_CASES}"
  run_suite "nstep${nsteps}_1rank_cpu" "${nsteps}" 1 "${ONE_RANK_CPU_CASES}"
  run_suite "nstep${nsteps}_${CPU_RANKS}rank_cpu_ref" "${nsteps}" "${CPU_RANKS}" \
    "${CPU_CASES}"
done

log "ALL DONE root=${FOLLOWUP_ROOT}"
[[ -f "${COMBINED}" ]] && log "combined=${COMBINED}"
