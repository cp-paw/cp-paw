#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

TEST=${TEST:-si64}
NSTEPS_LIST=${NSTEPS_LIST:-${NSTEPS:-2}}
EMPTY_BANDS_LIST=${EMPTY_BANDS_LIST:-${EMPTY_BANDS:-512}}
GPU_RANKS=${GPU_RANKS:-1}
SHARED_GPU_RANKS=${SHARED_GPU_RANKS:-4}
SHARED_EMPTY_BANDS=${SHARED_EMPTY_BANDS:-512}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-900}
DRY_RUN=${DRY_RUN:-no}
RUN_SHARED_GPU=${RUN_SHARED_GPU:-no}
RUN_LARGE_GPU=${RUN_LARGE_GPU:-no}
LARGE_EMPTY_BANDS=${LARGE_EMPTY_BANDS:-2048}
COPY_TOP=${COPY_TOP:-8}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-${COPY_TOP}}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-${COPY_TOP}}
RUN_ROOT_BASE=${RUN_ROOT_BASE:-${RUN_ROOT:-"${HERE}/runs/psim-lifecycle-$(date +%Y%m%d-%H%M%S)"}}
PSIM_LIFECYCLE_CASES=${PSIM_LIFECYCLE_CASES:-"gpu_resident gpu_psim_propagate gpu_resident_psim_phase gpu_resident_hpsi gpu_hpsi_psim_propagate gpu_resident_hpsi_psim_phase gpu_resident_hpsi_opsi gpu_resident_hpsi_opsi_psim_phase"}

COMBINED="${RUN_ROOT_BASE}-combined.tsv"
COPY_COMBINED="${RUN_ROOT_BASE}-copy-rows.md"
COPY_COMBINED_TSV="${RUN_ROOT_BASE}-copy-rows.tsv"
TRANSFER_COMBINED="${RUN_ROOT_BASE}-combined_transfer_rows.md"
TRANSFER_COMBINED_TSV="${RUN_ROOT_BASE}-combined_transfer_rows.tsv"
PRESENT_COMBINED="${RUN_ROOT_BASE}-combined_present_rows.md"
PRESENT_COMBINED_TSV="${RUN_ROOT_BASE}-combined_present_rows.tsv"
: > "${COMBINED}"

declare -a SUITE_ROOTS=()

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
  local label=$1
  local nsteps=$2
  local ranks=$3
  local empty_bands=$4
  local timeout=$5
  local root="${RUN_ROOT_BASE}-${label}"

  echo "Running PSIM lifecycle suite: ${label}"
  TEST="${TEST}" NSTEPS="${nsteps}" EMPTY_BANDS="${empty_bands}" \
    RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${PSIM_LIFECYCLE_CASES}" \
    TIMEOUT="${timeout}" RUN_ROOT="${root}" DRY_RUN="${DRY_RUN}" \
    REQUIRE_CASES="yes" \
    ENERGY_CHECK="${ENERGY_CHECK:-no}" \
    "${HERE}/run_benchmark.sh"
  append_suite "${label}" "${root}/benchmark.tsv"
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${COPY_TOP}" \
      --include-repeat --markdown "${root}" > "${root}/copy_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${COPY_TOP}" \
      --include-repeat "${root}" > "${root}/copy_rows.tsv"
  fi
  SUITE_ROOTS+=("${root}")
}

for nsteps in ${NSTEPS_LIST}; do
  for empty_bands in ${EMPTY_BANDS_LIST}; do
    run_suite "empty${empty_bands}-nstep${nsteps}-${GPU_RANKS}r" \
      "${nsteps}" "${GPU_RANKS}" "${empty_bands}" "${TIMEOUT}"
  done
done

case "${RUN_SHARED_GPU}" in
  yes|true|1)
    for nsteps in ${NSTEPS_LIST}; do
      run_suite "empty${SHARED_EMPTY_BANDS}-nstep${nsteps}-${SHARED_GPU_RANKS}r" \
        "${nsteps}" "${SHARED_GPU_RANKS}" "${SHARED_EMPTY_BANDS}" "${TIMEOUT}"
    done
    ;;
esac

case "${RUN_LARGE_GPU}" in
  yes|true|1)
    for nsteps in ${NSTEPS_LIST}; do
      run_suite "empty${LARGE_EMPTY_BANDS}-nstep${nsteps}-${GPU_RANKS}r" \
        "${nsteps}" "${GPU_RANKS}" "${LARGE_EMPTY_BANDS}" "${TIMEOUT}"
    done
    ;;
esac

if [[ -s "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${RUN_ROOT_BASE}-combined.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${RUN_ROOT_BASE}-compare.md" || true
  echo "Combined benchmark data: ${COMBINED}"
fi

if [[ "${#SUITE_ROOTS[@]}" -gt 0 ]]; then
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${COPY_TOP}" \
      --markdown "${SUITE_ROOTS[@]}" > "${COPY_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${COPY_TOP}" \
      "${SUITE_ROOTS[@]}" > "${COPY_COMBINED_TSV}" || true
    echo "Combined copy-row data: ${COPY_COMBINED}"
  else
    echo "Combined copy-row data: none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" --markdown \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" > "${TRANSFER_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" > "${TRANSFER_COMBINED_TSV}" || true
    echo "Combined transfer row data: ${TRANSFER_COMBINED}"
  else
    echo "Combined transfer row data: none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --markdown --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" > "${PRESENT_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" > "${PRESENT_COMBINED_TSV}" || true
    echo "Combined present row data: ${PRESENT_COMBINED}"
  else
    echo "Combined present row data: none"
  fi
fi
