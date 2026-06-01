#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
GPU_EXPLORATION_ROOT=${GPU_EXPLORATION_ROOT:-"${HERE}/runs/${TEST}-gpu-exploration-$(date +%Y%m%d-%H%M%S)"}
NSTEPS=${NSTEPS:-1}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
EMPTY_BANDS=${EMPTY_BANDS:-128}
DRY_RUN=${DRY_RUN:-no}
GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_stack gpu_resident_stack_cufft gpu_resident_nosync gpu gpu_force_all gpu_3dfft gpu_managed gpu_unified gpu_off"}
CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-16}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-16}

mkdir -p "${GPU_EXPLORATION_ROOT}"
echo "${GPU_EXPLORATION_ROOT}" > "${HERE}/runs/latest_gpu_exploration"

COMBINED="${GPU_EXPLORATION_ROOT}/combined_benchmark.tsv"
LOG="${GPU_EXPLORATION_ROOT}/gpu_exploration.log"
: > "${LOG}"

declare -a SUITE_ROOTS=()
FAILED_SUITES=0

log() {
  printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG}"
}

require_cases_for_suite() {
  case "${DRY_RUN}" in
    yes|true|1) echo yes ;;
    *) echo no ;;
  esac
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
  local require_cases

  if [[ -z ${cases// } ]]; then
    log "SKIP  suite=${suite} empty case list"
    echo "skipped" > "${root}.status"
    return 0
  fi

  require_cases=$(require_cases_for_suite)
  log "START suite=${suite} ranks=${ranks} cases=${cases}"
  if env TEST="${TEST}" EMPTY_BANDS="${EMPTY_BANDS}" NSTEPS="${NSTEPS}" \
      RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${cases}" TIMEOUT="${TIMEOUT}" \
      DRY_RUN="${DRY_RUN}" REQUIRE_CASES="${require_cases}" RUN_ROOT="${root}" \
      "${HERE}/run_benchmark.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${root}.status"
    FAILED_SUITES=$((FAILED_SUITES+1))
  fi
  append_suite "${suite}" "${root}/benchmark.tsv"
  SUITE_ROOTS+=("${root}")
}

if [[ -x "${HERE}/../../../src/Tools/Scripts/paw_gpu_capabilities.sh" ]]; then
  "${HERE}/../../../src/Tools/Scripts/paw_gpu_capabilities.sh" \
    > "${GPU_EXPLORATION_ROOT}/gpu_capabilities.txt" 2>&1 || true
fi

run_suite "one_rank_gpu" 1 "${GPU_CASES}"
run_suite "one_rank_cpu" 1 "${CPU_CASES}"
run_suite "eight_rank_cpu" 8 "${CPU_CASES}"

log "ALL DONE root=${GPU_EXPLORATION_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${GPU_EXPLORATION_ROOT}/combined_benchmark.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${GPU_EXPLORATION_ROOT}/combined_compare.md" || true
  log "combined=${COMBINED}"
fi

if [[ "${#SUITE_ROOTS[@]}" -gt 0 ]]; then
  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" --markdown \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${GPU_EXPLORATION_ROOT}/combined_transfer_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${GPU_EXPLORATION_ROOT}/combined_transfer_rows.tsv" || true
    log "transfer_rows=${GPU_EXPLORATION_ROOT}/combined_transfer_rows.md"
  else
    log "transfer_rows=none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --markdown --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${GPU_EXPLORATION_ROOT}/combined_present_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${GPU_EXPLORATION_ROOT}/combined_present_rows.tsv" || true
    log "present_rows=${GPU_EXPLORATION_ROOT}/combined_present_rows.md"
  else
    log "present_rows=none"
  fi
fi

if (( FAILED_SUITES > 0 )); then
  log "FAILED suites=${FAILED_SUITES}"
  case "${DRY_RUN}" in
    yes|true|1) exit 1 ;;
  esac
fi
