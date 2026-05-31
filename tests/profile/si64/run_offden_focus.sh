#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

TEST=${TEST:-si64}
NSTEPS=${NSTEPS:-1}
EMPTY_BANDS=${EMPTY_BANDS:-2048}
SHARED_EMPTY_BANDS=${SHARED_EMPTY_BANDS:-512}
GPU_RANKS=${GPU_RANKS:-1}
SHARED_GPU_RANKS=${SHARED_GPU_RANKS:-4}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-720}
RUN_SHARED_GPU=${RUN_SHARED_GPU:-yes}
RUN_ROOT_BASE=${RUN_ROOT_BASE:-${RUN_ROOT:-"${HERE}/runs/offden-focus-$(date +%Y%m%d-%H%M%S)"}}
OFFDEN_CASES=${OFFDEN_CASES:-"gpu_resident_hpsi_offden_cublas_devicepack gpu_resident_hpsi_offden_cublas_devicepack_accum gpu_resident_hpsi_offden_cublas_devicepack_proj gpu_resident_hpsi_offden_cublas_devicepack_proj_accum gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum"}

run_suite() {
  local label=$1
  local ranks=$2
  local empty_bands=$3
  local timeout=$4
  local root="${RUN_ROOT_BASE}-${label}"

  echo "Running off-site DENMAT focus suite: ${label}"
  TEST="${TEST}" NSTEPS="${NSTEPS}" EMPTY_BANDS="${empty_bands}" \
    RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${OFFDEN_CASES}" \
    TIMEOUT="${timeout}" RUN_ROOT="${root}" REQUIRE_CASES="yes" \
    "${HERE}/run_benchmark.sh"
}

run_suite "${EMPTY_BANDS}-${GPU_RANKS}r" "${GPU_RANKS}" "${EMPTY_BANDS}" "${TIMEOUT}"

case "${RUN_SHARED_GPU}" in
  yes|true|1)
    run_suite "${SHARED_EMPTY_BANDS}-${SHARED_GPU_RANKS}r" \
      "${SHARED_GPU_RANKS}" "${SHARED_EMPTY_BANDS}" "${TIMEOUT}"
    ;;
esac
