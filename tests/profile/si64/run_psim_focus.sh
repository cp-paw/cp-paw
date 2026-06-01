#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

TEST=${TEST:-si64}
NSTEPS=${NSTEPS:-1}
EMPTY_BANDS=${EMPTY_BANDS:-512}
SHARED_EMPTY_BANDS=${SHARED_EMPTY_BANDS:-512}
GPU_RANKS=${GPU_RANKS:-1}
SHARED_GPU_RANKS=${SHARED_GPU_RANKS:-4}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-720}
RUN_SHARED_GPU=${RUN_SHARED_GPU:-yes}
RUN_LARGE_GPU=${RUN_LARGE_GPU:-no}
LARGE_EMPTY_BANDS=${LARGE_EMPTY_BANDS:-2048}
RUN_ROOT_BASE=${RUN_ROOT_BASE:-${RUN_ROOT:-"${HERE}/runs/psim-focus-$(date +%Y%m%d-%H%M%S)"}}
PSIM_CASES=${PSIM_CASES:-"gpu_resident gpu_psim_propagate gpu_resident_psim_phase gpu_resident_hpsi gpu_hpsi_psim_propagate gpu_resident_hpsi_psim_phase"}

COMBINED="${RUN_ROOT_BASE}-combined.tsv"
: > "${COMBINED}"

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
  local ranks=$2
  local empty_bands=$3
  local timeout=$4
  local root="${RUN_ROOT_BASE}-${label}"

  echo "Running PSIM propagation focus suite: ${label}"
  TEST="${TEST}" NSTEPS="${NSTEPS}" EMPTY_BANDS="${empty_bands}" \
    RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${PSIM_CASES}" \
    TIMEOUT="${timeout}" RUN_ROOT="${root}" REQUIRE_CASES="yes" \
    "${HERE}/run_benchmark.sh"
  append_suite "${label}" "${root}/benchmark.tsv"
}

run_suite "${EMPTY_BANDS}-${GPU_RANKS}r" "${GPU_RANKS}" "${EMPTY_BANDS}" "${TIMEOUT}"

case "${RUN_SHARED_GPU}" in
  yes|true|1)
    run_suite "${SHARED_EMPTY_BANDS}-${SHARED_GPU_RANKS}r" \
      "${SHARED_GPU_RANKS}" "${SHARED_EMPTY_BANDS}" "${TIMEOUT}"
    ;;
esac

case "${RUN_LARGE_GPU}" in
  yes|true|1)
    run_suite "${LARGE_EMPTY_BANDS}-${GPU_RANKS}r" \
      "${GPU_RANKS}" "${LARGE_EMPTY_BANDS}" "${TIMEOUT}"
    ;;
esac

if [[ -s "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${RUN_ROOT_BASE}-combined.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${RUN_ROOT_BASE}-compare.md" || true
  echo "Combined benchmark data: ${COMBINED}"
fi
