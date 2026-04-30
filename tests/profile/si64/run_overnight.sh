#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
TEST=${TEST:-si64}
OVERNIGHT_ROOT=${OVERNIGHT_ROOT:-"${HERE}/runs/${TEST}-overnight-$(date +%Y%m%d-%H%M%S)"}
TIMEOUT=${TIMEOUT:-7200}
LONG_NSTEPS=${LONG_NSTEPS:-200}
SCALING_NSTEPS=${SCALING_NSTEPS:-100}
THRESHOLD_NSTEPS=${THRESHOLD_NSTEPS:-100}
NSYS_NSTEPS=${NSYS_NSTEPS:-5}
NVLAMATH_NSTEPS=${NVLAMATH_NSTEPS:-1}
CUFFTW_NSTEPS=${CUFFTW_NSTEPS:-1}
CUFFT_NSTEPS=${CUFFT_NSTEPS:-1}
GPU_ACC_NSTEPS=${GPU_ACC_NSTEPS:-1}
CUSOLVER_NSTEPS=${CUSOLVER_NSTEPS:-1}
BAND_NSTEPS=${BAND_NSTEPS:-3}
LONG_REPEATS=${LONG_REPEATS:-3}
SCALING_REPEATS=${SCALING_REPEATS:-2}
THRESHOLD_REPEATS=${THRESHOLD_REPEATS:-2}
BAND_REPEATS=${BAND_REPEATS:-1}
NSYS_RANKS=${NSYS_RANKS:-1}
RUN_NVLAMATH=${RUN_NVLAMATH:-no}
RUN_CUFFTW=${RUN_CUFFTW:-no}
RUN_CUFFT=${RUN_CUFFT:-no}
RUN_GPU_ACC=${RUN_GPU_ACC:-no}
RUN_GPU_DIAGNOSTICS=${RUN_GPU_DIAGNOSTICS:-no}
RUN_CUSOLVER=${RUN_CUSOLVER:-no}
RUN_BAND_BENCHMARK=${RUN_BAND_BENCHMARK:-no}
MAIN_CASES=${MAIN_CASES:-"nvpl cublas cublas_off"}
SCALING_CASES=${SCALING_CASES:-"nvpl cublas"}
GPU_ACC_CASES=${GPU_ACC_CASES:-"cpu nvpl gpu gpu_off"}
GPU_DIAGNOSTIC_CASES=${GPU_DIAGNOSTIC_CASES:-"gpu_force_all gpu_no_cufft gpu_no_cublas gpu_no_cusolver"}
BAND_TEST=${BAND_TEST:-si64_bands}
BAND_RANKS=${BAND_RANKS:-1}
BAND_EMPTY_BANDS=${BAND_EMPTY_BANDS:-128}
BAND_CASES=${BAND_CASES:-"cpu nvpl gpu gpu_off"}
THRESHOLDS=${THRESHOLDS:-"1e7"}

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}

nvhpc_platform() {
  case "$(uname -s)_$(uname -m)" in
    Linux_aarch64|Linux_arm64) echo "Linux_aarch64" ;;
    Linux_x86_64) echo "Linux_x86_64" ;;
    *) echo "" ;;
  esac
}

default_mpirun() {
  local platform root candidate
  platform=$(nvhpc_platform)
  for root in "${NVHPC_ROOT:-}" /opt/nvidia/hpc_sdk/${platform}/*; do
    [[ -n "${root}" && -d "${root}" ]] || continue
    for candidate in \
        "${root}/comm_libs/hpcx/bin/mpirun" \
        "${root}"/comm_libs/*/hpcx/*/ompi/bin/mpirun \
        "${root}"/comm_libs/*/hpcx/bin/mpirun; do
      if [[ -x "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    done
  done
  command -v mpirun 2>/dev/null || echo mpirun
}

export MPIRUN=${MPIRUN:-$(default_mpirun)}

mkdir -p "${OVERNIGHT_ROOT}"
echo "${OVERNIGHT_ROOT}" > "${HERE}/runs/latest_overnight"

COMBINED="${OVERNIGHT_ROOT}/combined_benchmark.tsv"
SUMMARY_LOG="${OVERNIGHT_ROOT}/overnight.log"
: > "${SUMMARY_LOG}"

log() {
  printf '%s %s\n' "$(date -Is)" "$*" | tee -a "${SUMMARY_LOG}"
}

label_safe() {
  printf '%s' "$1" | tr '.+' 'pp' | tr -c 'A-Za-z0-9_-' '_'
}

capture_metadata() {
  {
    echo "date=$(date -Is)"
    echo "hostname=$(hostname)"
    echo "root=${ROOT}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
    echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
    echo
    uname -a
    echo
    command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi || true
    echo
    command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi topo -m || true
    echo
    command -v nvaccelinfo >/dev/null 2>&1 && nvaccelinfo || true
    echo
    command -v nsys >/dev/null 2>&1 && nsys --version || true
    echo
    "${MPIRUN}" --version || true
    echo
    ls -l "${ROOT}/bin/nvhpc_profile/paw_nvhpc_profile.x" \
          "${ROOT}/bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" \
          "${ROOT}/bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" \
          "${ROOT}/bin/nvhpc_nvlamath_profile_parallel/ppaw_nvhpc_nvlamath_profile.x" \
          "${ROOT}/bin/nvhpc_cufft_profile/paw_nvhpc_cufft_profile.x" \
          "${ROOT}/bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" \
          "${ROOT}/bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x" \
          "${ROOT}/bin/nvhpc_gpu_acc_profile_parallel/ppaw_nvhpc_gpu_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cufft_cublas_acc_profile/paw_nvhpc_cufft_cublas_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cufft_cublas_acc_profile_parallel/ppaw_nvhpc_cufft_cublas_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cublas_acc_profile/paw_nvhpc_cublas_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cusolver_acc_profile/paw_nvhpc_cusolver_acc_profile.x" \
          "${ROOT}/bin/nvhpc_cusolver_acc_profile_parallel/ppaw_nvhpc_cusolver_acc_profile.x" || true
  } > "${OVERNIGHT_ROOT}/metadata.txt" 2>&1
}

append_suite() {
  local suite=$1
  local tsv=$2
  [[ -f "${tsv}" ]] || return 0
  if [[ ! -s "${COMBINED}" ]]; then
    awk 'NR == 1 { print "suite\t" $0; next } { print suite "\t" $0 }' suite="${suite}" "${tsv}" > "${COMBINED}"
  else
    awk 'NR > 1 { print suite "\t" $0 }' suite="${suite}" "${tsv}" >> "${COMBINED}"
  fi
}

run_suite() {
  local suite=$1
  local nsteps=$2
  local ranks=$3
  local repeats=$4
  local cases=$5
  shift 5
  local suite_root="${OVERNIGHT_ROOT}/${suite}"
  local suite_log="${OVERNIGHT_ROOT}/${suite}.log"

  log "START suite=${suite} nsteps=${nsteps} ranks=${ranks} repeats=${repeats} cases=${cases} env=[$*]"
  if env NSTEPS="${nsteps}" RANKS="${ranks}" REPEATS="${repeats}" CASES="${cases}" \
      TIMEOUT="${TIMEOUT}" RUN_ROOT="${suite_root}" "$@" \
      "${HERE}/run_benchmark.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${suite_root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${suite_root}.status"
  fi
  append_suite "${suite}" "${suite_root}/benchmark.tsv"
}

run_nsys_trace() {
  local suite="nsys_nstep${NSYS_NSTEPS}_${NSYS_RANKS}ranks"
  local suite_root="${OVERNIGHT_ROOT}/${suite}"
  local suite_log="${OVERNIGHT_ROOT}/${suite}.log"

  log "START suite=${suite}"
  if env NSTEPS="${NSYS_NSTEPS}" RANKS="${NSYS_RANKS}" TIMEOUT="${TIMEOUT}" \
      RUN_ROOT="${suite_root}" "${HERE}/run_nsys.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${suite_root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${suite_root}.status"
  fi
}

capture_metadata

run_suite "main_${LONG_NSTEPS}steps_4ranks" "${LONG_NSTEPS}" 4 "${LONG_REPEATS}" \
  "${MAIN_CASES}"

for ranks in 1 2 4; do
  run_suite "scaling_${SCALING_NSTEPS}steps_${ranks}ranks" "${SCALING_NSTEPS}" "${ranks}" \
    "${SCALING_REPEATS}" "${SCALING_CASES}"
done

case "${RUN_NVLAMATH}" in
  yes|true|1)
    run_suite "nvlamath_${NVLAMATH_NSTEPS}steps_1rank" "${NVLAMATH_NSTEPS}" 1 1 \
      "cpu nvpl nvlamath"
    ;;
esac

case "${RUN_CUFFTW}" in
  yes|true|1)
    run_suite "cufftw_${CUFFTW_NSTEPS}steps_1rank" "${CUFFTW_NSTEPS}" 1 1 \
      "nvpl cufftw"
    ;;
esac

case "${RUN_CUFFT}" in
  yes|true|1)
    run_suite "cufft_${CUFFT_NSTEPS}steps_4ranks" "${CUFFT_NSTEPS}" 4 1 \
      "nvpl cufft cufft_off"
    ;;
esac

case "${RUN_GPU_ACC}" in
  yes|true|1)
    gpu_acc_cases="${GPU_ACC_CASES}"
    case "${RUN_GPU_DIAGNOSTICS}" in
      yes|true|1) gpu_acc_cases="${gpu_acc_cases} ${GPU_DIAGNOSTIC_CASES}" ;;
    esac
    run_suite "gpu_acc_${GPU_ACC_NSTEPS}steps_1rank" "${GPU_ACC_NSTEPS}" 1 1 \
      "${gpu_acc_cases}"
    run_suite "cpu_ref_${GPU_ACC_NSTEPS}steps_8ranks" "${GPU_ACC_NSTEPS}" 8 1 \
      "cpu nvpl"
    ;;
esac

case "${RUN_CUSOLVER}" in
  yes|true|1)
    run_suite "cusolver_${CUSOLVER_NSTEPS}steps_1rank" "${CUSOLVER_NSTEPS}" 1 1 \
      "cpu nvpl cusolver cusolver_conservative cusolver_off"
    run_suite "cusolver_cpu_ref_${CUSOLVER_NSTEPS}steps_8ranks" "${CUSOLVER_NSTEPS}" 8 1 \
      "cpu nvpl"
    ;;
esac

case "${RUN_BAND_BENCHMARK}" in
  yes|true|1)
    run_suite "${BAND_TEST}_${BAND_NSTEPS}steps_${BAND_RANKS}ranks" \
      "${BAND_NSTEPS}" "${BAND_RANKS}" "${BAND_REPEATS}" "${BAND_CASES}" \
      "TEST=${BAND_TEST}" "EMPTY_BANDS=${BAND_EMPTY_BANDS}"
    ;;
esac

for threshold in ${THRESHOLDS}; do
  safe_threshold=$(label_safe "${threshold}")
  run_suite "threshold_${safe_threshold}_${THRESHOLD_NSTEPS}steps_4ranks" "${THRESHOLD_NSTEPS}" 4 \
    "${THRESHOLD_REPEATS}" "cublas" "CPPAW_CUBLAS_ACC_MINFLOP=${threshold}"
done

run_nsys_trace

log "ALL DONE root=${OVERNIGHT_ROOT}"
[[ -f "${COMBINED}" ]] && log "combined=${COMBINED}"
