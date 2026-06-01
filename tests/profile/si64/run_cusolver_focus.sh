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
DRY_RUN=${DRY_RUN:-no}
CUSOLVER_CASES=${CUSOLVER_CASES:-"cusolver cusolver_generalized cusolver_generalized_conservative cusolver_off"}
CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}
SOLVER_ROW_TOP=${SOLVER_ROW_TOP:-16}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-16}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-16}

mkdir -p "${CUSOLVER_FOCUS_ROOT}"
echo "${CUSOLVER_FOCUS_ROOT}" > "${HERE}/runs/latest_cusolver_focus"

COMBINED="${CUSOLVER_FOCUS_ROOT}/combined_benchmark.tsv"
LOG="${CUSOLVER_FOCUS_ROOT}/cusolver_focus.log"
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
  local empty_bands=$2
  local ranks=$3
  local cases=$4
  local root="${CUSOLVER_FOCUS_ROOT}/${suite}"
  local suite_log="${CUSOLVER_FOCUS_ROOT}/${suite}.log"
  local require_cases

  if [[ -z ${cases// } ]]; then
    log "SKIP  suite=${suite} empty case list"
    echo "skipped" > "${root}.status"
    return 0
  fi

  require_cases=$(require_cases_for_suite)
  log "START suite=${suite} empty_bands=${empty_bands} ranks=${ranks} cases=${cases}"
  if env TEST="${TEST}" EMPTY_BANDS="${empty_bands}" NSTEPS="${NSTEPS}" \
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

if [[ "${#SUITE_ROOTS[@]}" -gt 0 ]]; then
  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${SOLVER_ROW_TOP}" --markdown --include-zero --sort-by seconds \
      --op-prefix LAPACK --op-prefix CUSOLVER \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_solver_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${SOLVER_ROW_TOP}" --include-zero --sort-by seconds \
      --op-prefix LAPACK --op-prefix CUSOLVER \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_solver_rows.tsv" || true
    log "solver_rows=${CUSOLVER_FOCUS_ROOT}/combined_solver_rows.md"
  else
    log "solver_rows=none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" --markdown \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_transfer_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_transfer_rows.tsv" || true
    log "transfer_rows=${CUSOLVER_FOCUS_ROOT}/combined_transfer_rows.md"
  else
    log "transfer_rows=none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --markdown --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_present_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${CUSOLVER_FOCUS_ROOT}/combined_present_rows.tsv" || true
    log "present_rows=${CUSOLVER_FOCUS_ROOT}/combined_present_rows.md"
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
