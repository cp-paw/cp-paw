#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

TEST=${TEST:-si64_bands}
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
ROW_TOP=${ROW_TOP:-16}
VPSI_ROW_TOP=${VPSI_ROW_TOP:-8}
FFT_ROW_TOP=${FFT_ROW_TOP:-12}
RUN_ROOT_BASE=${RUN_ROOT_BASE:-${RUN_ROOT:-"${HERE}/runs/vpsi-boundary-$(date +%Y%m%d-%H%M%S)"}}
VPSI_BOUNDARY_CASES=${VPSI_BOUNDARY_CASES:-"gpu_resident_hpsi gpu_resident_hpsi_opsi gpu_resident_stack gpu_resident_stack_cufft gpu_resident_stack_cufft_force gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_accmap"}

COMBINED="${RUN_ROOT_BASE}-combined.tsv"
ROWS_COMBINED="${RUN_ROOT_BASE}-profile-rows.md"
ROWS_COMBINED_TSV="${RUN_ROOT_BASE}-profile-rows.tsv"
VPSI_ROWS_COMBINED="${RUN_ROOT_BASE}-vpsi-rows.md"
VPSI_ROWS_COMBINED_TSV="${RUN_ROOT_BASE}-vpsi-rows.tsv"
FFT_ROWS_COMBINED="${RUN_ROOT_BASE}-fft-phase-rows.md"
FFT_ROWS_COMBINED_TSV="${RUN_ROOT_BASE}-fft-phase-rows.tsv"
: > "${COMBINED}"

declare -a SUITE_ROOTS=()
PROFILE_ROW_ARGS=(
  --op-prefix ACC_COPY
  --op-prefix ACC_UPDATE
  --op-prefix ACC_PRESENT
  --include-zero
)
VPSI_ROW_ARGS=(
  --op-prefix PAW_VPSI_
  --include-zero
  --sort-by seconds
)
FFT_ROW_ARGS=(
  --op-prefix PW_FFT_
  --op-prefix PW_GTOR_
  --op-prefix PW_RTOG_
  --include-zero
  --sort-by seconds
)

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

  echo "Running VPSI boundary suite: ${label}"
  TEST="${TEST}" NSTEPS="${nsteps}" EMPTY_BANDS="${empty_bands}" \
    RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${VPSI_BOUNDARY_CASES}" \
    TIMEOUT="${timeout}" RUN_ROOT="${root}" DRY_RUN="${DRY_RUN}" \
    REQUIRE_CASES="yes" \
    ENERGY_CHECK="${ENERGY_CHECK:-no}" \
    "${HERE}/run_benchmark.sh"
  append_suite "${label}" "${root}/benchmark.tsv"
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${ROW_TOP}" \
      --include-repeat --markdown "${PROFILE_ROW_ARGS[@]}" "${root}" \
      > "${root}/profile_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${ROW_TOP}" \
      --include-repeat "${PROFILE_ROW_ARGS[@]}" "${root}" \
      > "${root}/profile_rows.tsv"
  fi
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${VPSI_ROW_TOP}" \
      --include-repeat --markdown "${VPSI_ROW_ARGS[@]}" "${root}" \
      > "${root}/vpsi_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${VPSI_ROW_TOP}" \
      --include-repeat "${VPSI_ROW_ARGS[@]}" "${root}" \
      > "${root}/vpsi_rows.tsv"
  fi
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${FFT_ROW_TOP}" \
      --include-repeat --markdown "${FFT_ROW_ARGS[@]}" "${root}" \
      > "${root}/fft_phase_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${FFT_ROW_TOP}" \
      --include-repeat "${FFT_ROW_ARGS[@]}" "${root}" \
      > "${root}/fft_phase_rows.tsv"
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
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${ROW_TOP}" \
      --markdown "${PROFILE_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${ROWS_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${ROW_TOP}" \
      "${PROFILE_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${ROWS_COMBINED_TSV}" || true
    echo "Combined profile-row data: ${ROWS_COMBINED}"
  else
    echo "Combined profile-row data: none"
  fi
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${VPSI_ROW_TOP}" \
      --markdown "${VPSI_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${VPSI_ROWS_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${VPSI_ROW_TOP}" \
      "${VPSI_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${VPSI_ROWS_COMBINED_TSV}" || true
    echo "Combined VPSI-row data: ${VPSI_ROWS_COMBINED}"
  else
    echo "Combined VPSI-row data: none"
  fi
  if python3 "${HERE}/profile_copy_rows.py" --per-case --top "${FFT_ROW_TOP}" \
      --markdown "${FFT_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${FFT_ROWS_COMBINED}"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case --top "${FFT_ROW_TOP}" \
      "${FFT_ROW_ARGS[@]}" "${SUITE_ROOTS[@]}" \
      > "${FFT_ROWS_COMBINED_TSV}" || true
    echo "Combined FFT phase-row data: ${FFT_ROWS_COMBINED}"
  else
    echo "Combined FFT phase-row data: none"
  fi
fi
