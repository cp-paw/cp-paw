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
EMPTY_BANDS_LIST=${EMPTY_BANDS_LIST:-${EMPTY_BANDS}}
DRY_RUN=${DRY_RUN:-no}
GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_stack gpu_resident_nosync gpu gpu_off"}
CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}
ONE_RANK_CPU_CASES=${ONE_RANK_CPU_CASES-"cpu nvhpc_cpu"}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-16}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-16}

mkdir -p "${FOLLOWUP_ROOT}"
echo "${FOLLOWUP_ROOT}" > "${HERE}/runs/latest_followup"

COMBINED="${FOLLOWUP_ROOT}/combined_benchmark.tsv"
SUMMARY_LOG="${FOLLOWUP_ROOT}/followup.log"
: > "${SUMMARY_LOG}"

declare -a SUITE_ROOTS=()
FAILED_SUITES=0

log() {
  printf '%s %s\n' "$(iso_now)" "$*" | tee -a "${SUMMARY_LOG}"
}

iso_now() {
  date -u +%Y-%m-%dT%H:%M:%SZ
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
  local nsteps=$2
  local ranks=$3
  local cases=$4
  shift 4
  local suite_root="${FOLLOWUP_ROOT}/${suite}"
  local suite_log="${FOLLOWUP_ROOT}/${suite}.log"
  local require_cases

  if [[ -z ${cases// } ]]; then
    log "SKIP  suite=${suite} empty case list"
    echo "skipped" > "${suite_root}.status"
    return 0
  fi

  require_cases=$(require_cases_for_suite)
  log "START suite=${suite} nsteps=${nsteps} ranks=${ranks} cases=${cases} env=[$*]"
  if env TEST="${TEST}" EMPTY_BANDS="${EMPTY_BANDS}" NSTEPS="${nsteps}" \
      RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${cases}" TIMEOUT="${TIMEOUT}" \
      DRY_RUN="${DRY_RUN}" REQUIRE_CASES="${require_cases}" RUN_ROOT="${suite_root}" \
      "$@" "${HERE}/run_benchmark.sh" \
      > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${suite_root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${suite_root}.status"
    FAILED_SUITES=$((FAILED_SUITES+1))
  fi
  append_suite "${suite}" "${suite_root}/benchmark.tsv"
  SUITE_ROOTS+=("${suite_root}")
}

for empty_bands in ${EMPTY_BANDS_LIST}; do
  for nsteps in ${NSTEPS_LIST}; do
    run_suite "empty${empty_bands}_nstep${nsteps}_${GPU_RANKS}rank_gpu" \
      "${nsteps}" "${GPU_RANKS}" "${GPU_CASES}" "EMPTY_BANDS=${empty_bands}"
    run_suite "empty${empty_bands}_nstep${nsteps}_1rank_cpu" \
      "${nsteps}" 1 "${ONE_RANK_CPU_CASES}" "EMPTY_BANDS=${empty_bands}"
    run_suite "empty${empty_bands}_nstep${nsteps}_${CPU_RANKS}rank_cpu_ref" \
      "${nsteps}" "${CPU_RANKS}" "${CPU_CASES}" "EMPTY_BANDS=${empty_bands}"
  done
done

log "ALL DONE root=${FOLLOWUP_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${FOLLOWUP_ROOT}/combined_benchmark.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${FOLLOWUP_ROOT}/combined_compare.md" || true
  log "combined=${COMBINED}"
fi

if [[ "${#SUITE_ROOTS[@]}" -gt 0 ]]; then
  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" --markdown \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${FOLLOWUP_ROOT}/combined_transfer_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${FOLLOWUP_ROOT}/combined_transfer_rows.tsv" || true
    log "transfer_rows=${FOLLOWUP_ROOT}/combined_transfer_rows.md"
  else
    log "transfer_rows=none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --markdown --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${FOLLOWUP_ROOT}/combined_present_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${FOLLOWUP_ROOT}/combined_present_rows.tsv" || true
    log "present_rows=${FOLLOWUP_ROOT}/combined_present_rows.md"
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
