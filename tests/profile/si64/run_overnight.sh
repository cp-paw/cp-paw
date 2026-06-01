#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
. "${HERE}/case_recommendations.sh"
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
BAND_NSTEPS_LIST=${BAND_NSTEPS_LIST:-${BAND_NSTEPS}}
LONG_REPEATS=${LONG_REPEATS:-3}
SCALING_REPEATS=${SCALING_REPEATS:-2}
THRESHOLD_REPEATS=${THRESHOLD_REPEATS:-2}
BAND_REPEATS=${BAND_REPEATS:-1}
NSYS_RANKS=${NSYS_RANKS:-1}
DRY_RUN=${DRY_RUN:-no}
RUN_NVLAMATH=${RUN_NVLAMATH:-no}
RUN_CUFFTW=${RUN_CUFFTW:-no}
RUN_CUFFT=${RUN_CUFFT:-no}
RUN_GPU_ACC=${RUN_GPU_ACC:-no}
RUN_GPU_DIAGNOSTICS=${RUN_GPU_DIAGNOSTICS:-no}
RUN_CUSOLVER=${RUN_CUSOLVER:-no}
RUN_BAND_BENCHMARK=${RUN_BAND_BENCHMARK:-no}
RUN_NSYS=${RUN_NSYS:-auto}
MAIN_CASES=${MAIN_CASES:-auto}
SCALING_CASES=${SCALING_CASES:-auto}
GPU_ACC_CASES=${GPU_ACC_CASES:-auto}
GPU_DIAGNOSTIC_CASES=${GPU_DIAGNOSTIC_CASES:-auto}
BAND_TEST=${BAND_TEST:-si64_bands}
BAND_RANKS=${BAND_RANKS:-1}
BAND_CPU_RANKS=${BAND_CPU_RANKS:-8}
BAND_EMPTY_BANDS=${BAND_EMPTY_BANDS:-128}
BAND_EMPTY_BANDS_LIST=${BAND_EMPTY_BANDS_LIST:-${BAND_EMPTY_BANDS}}
BAND_CASES=${BAND_CASES:-auto}
BAND_ONE_RANK_CPU_CASES=${BAND_ONE_RANK_CPU_CASES:-auto}
BAND_CPU_CASES=${BAND_CPU_CASES:-auto}
THRESHOLDS=${THRESHOLDS:-"1e7"}
THRESHOLD_CASES=${THRESHOLD_CASES:-auto}
NSYS_CASE=${NSYS_CASE:-auto}
PROFILE_ROW_TOP=${PROFILE_ROW_TOP:-16}
PRESENT_ROW_TOP=${PRESENT_ROW_TOP:-16}

DEFAULT_MAIN_CASES="nvhpc_cpu cublas cublas_off"
DEFAULT_SCALING_CASES="nvhpc_cpu cublas"
DEFAULT_GPU_ACC_CASES="cpu nvhpc_cpu gpu_resident_stack gpu_resident_off"
DEFAULT_GPU_DIAGNOSTIC_CASES="gpu_nosync gpu_resident_stack_cufft gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_accmap gpu_force_all gpu_3dfft gpu_no_cufft gpu_no_cublas gpu_no_cusolver"
DEFAULT_BAND_CASES="gpu_resident_stack gpu_resident_off"
DEFAULT_CPU_CASES="cpu nvhpc_cpu"
DEFAULT_THRESHOLD_CASES="cublas"
DEFAULT_NSYS_CASE="gpu_resident_stack"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export BLIS_NUM_THREADS=${BLIS_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NVPL_NUM_THREADS=${NVPL_NUM_THREADS:-1}

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
  for root in \
      "${NVHPC_ROOT:-}" \
      /opt/nvidia/hpc_sdk/${platform}/* \
      "${HOME:-}"/opt/nvidia/hpc_sdk/${platform}/*; do
    [[ -n "${root}" && -d "${root}" ]] || continue
    for candidate in \
        "${root}/comm_libs/hpcx/bin/mpirun" \
        "${root}"/comm_libs/*/hpcx/*/ompi/bin/mpirun \
        "${root}"/comm_libs/*/hpcx/bin/mpirun \
        "${root}"/../*/comm_libs/hpcx/bin/mpirun \
        "${root}"/../*/comm_libs/*/hpcx/*/ompi/bin/mpirun \
        "${root}"/../*/comm_libs/*/hpcx/bin/mpirun; do
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

declare -a SUITE_ROOTS=()

log() {
  printf '%s %s\n' "$(iso_now)" "$*" | tee -a "${SUMMARY_LOG}"
}

iso_now() {
  date -u +%Y-%m-%dT%H:%M:%SZ
}

first_case() {
  local item
  for item in "$@"; do
    [[ -n "${item}" ]] || continue
    echo "${item}"
    return 0
  done
}

label_safe() {
  printf '%s' "$1" | tr '.+' 'pp' | tr -c 'A-Za-z0-9_-' '_'
}

resolve_case_list() {
  local selected=$1
  local key=$2
  local fallback=$3
  local cases

  cases=$(cppaw_resolve_recommended_cases "${selected}" "${key}" "${fallback}")
  cppaw_dedup_space_list ${cases}
}

resolve_first_case() {
  local selected=$1
  local key=$2
  local fallback=$3
  local cases

  cases=$(resolve_case_list "${selected}" "${key}" "${fallback}")
  first_case ${cases}
}

resolve_threshold_cases() {
  local selected=$1
  local key=$2
  local fallback=$3

  case "${selected}" in
    auto|recommended)
      resolve_first_case "${selected}" "${key}" "${fallback}"
      ;;
    *)
      resolve_case_list "${selected}" "${key}" "${fallback}"
      ;;
  esac
}

run_enabled() {
  case "$1" in
    yes|true|1) return 0 ;;
    no|false|0) return 1 ;;
    auto|recommended) [[ -n ${2// } ]] ;;
    *) return 1 ;;
  esac
}

capture_metadata() {
  {
    echo "date=$(iso_now)"
    echo "hostname=$(hostname)"
    echo "root=${ROOT}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
    echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
    echo "BLIS_NUM_THREADS=${BLIS_NUM_THREADS}"
    echo "VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS}"
    echo "NVPL_NUM_THREADS=${NVPL_NUM_THREADS}"
    echo "DRY_RUN=${DRY_RUN}"
    echo "MAIN_CASES=${MAIN_CASES}"
    echo "SCALING_CASES=${SCALING_CASES}"
    echo "GPU_ACC_CASES=${GPU_ACC_CASES}"
    echo "GPU_DIAGNOSTIC_CASES=${GPU_DIAGNOSTIC_CASES}"
    echo "BAND_CASES=${BAND_CASES}"
    echo "BAND_ONE_RANK_CPU_CASES=${BAND_ONE_RANK_CPU_CASES}"
    echo "BAND_CPU_CASES=${BAND_CPU_CASES}"
    echo "THRESHOLD_CASES=${THRESHOLD_CASES}"
    echo "NSYS_CASE=${NSYS_CASE}"
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
    if [[ -x "${ROOT}/src/Tools/Scripts/paw_gpu_capabilities.sh" ]]; then
      "${ROOT}/src/Tools/Scripts/paw_gpu_capabilities.sh" || true
      echo
    fi
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
          "${ROOT}/bin/nvhpc_gpu_all_profile/paw_nvhpc_gpu_all_profile.x" \
          "${ROOT}/bin/nvhpc_gpu_all_profile_parallel/ppaw_nvhpc_gpu_all_profile.x" \
          "${ROOT}/bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" \
          "${ROOT}/bin/nvhpc_gpu_acc_residency_profile_parallel/ppaw_nvhpc_gpu_acc_residency_profile.x" \
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

  if [[ -z ${cases// } ]]; then
    log "SKIP  suite=${suite} empty case list"
    echo "skipped" > "${suite_root}.status"
    return 0
  fi

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
  SUITE_ROOTS+=("${suite_root}")
}

run_nsys_trace() {
  local nsys_case=$1
  local suite="nsys_nstep${NSYS_NSTEPS}_${NSYS_RANKS}ranks"
  local suite_root="${OVERNIGHT_ROOT}/${suite}"
  local suite_log="${OVERNIGHT_ROOT}/${suite}.log"

  if [[ -z ${nsys_case// } ]]; then
    log "SKIP  suite=${suite} empty Nsight case"
    echo "skipped" > "${suite_root}.status"
    return 0
  fi
  case "${DRY_RUN}" in
    yes|true|1)
      log "DRY-RUN skip suite=${suite} case=${nsys_case}"
      echo "dry-run" > "${suite_root}.status"
      return 0
      ;;
  esac

  log "START suite=${suite} case=${nsys_case}"
  if env NSTEPS="${NSYS_NSTEPS}" RANKS="${NSYS_RANKS}" TIMEOUT="${TIMEOUT}" \
      CASE="${nsys_case}" RUN_ROOT="${suite_root}" \
      "${HERE}/run_nsys.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${suite_root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${suite_root}.status"
  fi
}

if cppaw_write_capabilities_file "${OVERNIGHT_ROOT}/gpu_capabilities.txt"; then
  log "gpu_capabilities=${OVERNIGHT_ROOT}/gpu_capabilities.txt"
fi

MAIN_CASES=$(resolve_case_list "${MAIN_CASES}" recommended_resource_cases "${DEFAULT_MAIN_CASES}")
SCALING_CASES=$(resolve_case_list "${SCALING_CASES}" recommended_resource_cases "${DEFAULT_SCALING_CASES}")
GPU_ACC_CASES=$(resolve_case_list "${GPU_ACC_CASES}" recommended_resource_cases "${DEFAULT_GPU_ACC_CASES}")
GPU_DIAGNOSTIC_CASES=$(resolve_case_list "${GPU_DIAGNOSTIC_CASES}" recommended_gpu_diagnostic_cases "${DEFAULT_GPU_DIAGNOSTIC_CASES}")
BAND_CASES=$(resolve_case_list "${BAND_CASES}" recommended_gpu_cases "${DEFAULT_BAND_CASES}")
BAND_ONE_RANK_CPU_CASES=$(resolve_case_list "${BAND_ONE_RANK_CPU_CASES}" recommended_cpu_cases "${DEFAULT_CPU_CASES}")
BAND_CPU_CASES=$(resolve_case_list "${BAND_CPU_CASES}" recommended_cpu_cases "${DEFAULT_CPU_CASES}")
THRESHOLD_CASES=$(resolve_threshold_cases "${THRESHOLD_CASES}" recommended_gpu_cases "${DEFAULT_THRESHOLD_CASES}")
NSYS_CASE=$(resolve_first_case "${NSYS_CASE}" recommended_gpu_cases "${DEFAULT_NSYS_CASE}")

log "selected_cases main='${MAIN_CASES:-none}' scaling='${SCALING_CASES:-none}' gpu_acc='${GPU_ACC_CASES:-none}' diagnostics='${GPU_DIAGNOSTIC_CASES:-none}' band_gpu='${BAND_CASES:-none}' band_cpu='${BAND_CPU_CASES:-none}' threshold='${THRESHOLD_CASES:-none}' nsys='${NSYS_CASE:-none}'"

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
      "cpu nvhpc_cpu nvlamath"
    ;;
esac

case "${RUN_CUFFTW}" in
  yes|true|1)
    run_suite "cufftw_${CUFFTW_NSTEPS}steps_1rank" "${CUFFTW_NSTEPS}" 1 1 \
      "nvhpc_cpu cufftw"
    ;;
esac

case "${RUN_CUFFT}" in
  yes|true|1)
    run_suite "cufft_${CUFFT_NSTEPS}steps_4ranks" "${CUFFT_NSTEPS}" 4 1 \
      "nvhpc_cpu cufft cufft_off"
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
      "cpu nvhpc_cpu"
    ;;
esac

case "${RUN_CUSOLVER}" in
  yes|true|1)
    run_suite "cusolver_${CUSOLVER_NSTEPS}steps_1rank" "${CUSOLVER_NSTEPS}" 1 1 \
      "cpu nvhpc_cpu cusolver cusolver_conservative cusolver_off"
    run_suite "cusolver_cpu_ref_${CUSOLVER_NSTEPS}steps_8ranks" "${CUSOLVER_NSTEPS}" 8 1 \
      "cpu nvhpc_cpu"
    ;;
esac

case "${RUN_BAND_BENCHMARK}" in
  yes|true|1)
    for band_empty_bands in ${BAND_EMPTY_BANDS_LIST}; do
      for band_nsteps in ${BAND_NSTEPS_LIST}; do
        run_suite "${BAND_TEST}_empty${band_empty_bands}_${band_nsteps}steps_${BAND_RANKS}ranks_gpu" \
          "${band_nsteps}" "${BAND_RANKS}" "${BAND_REPEATS}" "${BAND_CASES}" \
          "TEST=${BAND_TEST}" "EMPTY_BANDS=${band_empty_bands}"
        run_suite "${BAND_TEST}_empty${band_empty_bands}_${band_nsteps}steps_1rank_cpu" \
          "${band_nsteps}" 1 "${BAND_REPEATS}" "${BAND_ONE_RANK_CPU_CASES}" \
          "TEST=${BAND_TEST}" "EMPTY_BANDS=${band_empty_bands}"
        run_suite "${BAND_TEST}_empty${band_empty_bands}_${band_nsteps}steps_${BAND_CPU_RANKS}ranks_cpu_ref" \
          "${band_nsteps}" "${BAND_CPU_RANKS}" "${BAND_REPEATS}" "${BAND_CPU_CASES}" \
          "TEST=${BAND_TEST}" "EMPTY_BANDS=${band_empty_bands}"
      done
    done
    ;;
esac

if [[ -n ${THRESHOLD_CASES// } ]]; then
  for threshold in ${THRESHOLDS}; do
    safe_threshold=$(label_safe "${threshold}")
    run_suite "threshold_${safe_threshold}_${THRESHOLD_NSTEPS}steps_4ranks" "${THRESHOLD_NSTEPS}" 4 \
      "${THRESHOLD_REPEATS}" "${THRESHOLD_CASES}" "CPPAW_CUBLAS_ACC_MINFLOP=${threshold}"
  done
else
  log "SKIP  suite=threshold empty case list"
fi

if run_enabled "${RUN_NSYS}" "${NSYS_CASE}"; then
  run_nsys_trace "${NSYS_CASE}"
else
  log "SKIP  suite=nsys disabled_or_empty_case"
fi

log "ALL DONE root=${OVERNIGHT_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${OVERNIGHT_ROOT}/combined_benchmark.md" || true
  python3 "${HERE}/benchmark_compare.py" "${COMBINED}" \
    > "${OVERNIGHT_ROOT}/combined_compare.md" || true
  log "combined=${COMBINED}"
fi

if [[ "${#SUITE_ROOTS[@]}" -gt 0 ]]; then
  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" --markdown \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${OVERNIGHT_ROOT}/combined_transfer_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PROFILE_ROW_TOP}" \
      --op-prefix ACC_COPY --op-prefix ACC_UPDATE \
      "${SUITE_ROOTS[@]}" \
      > "${OVERNIGHT_ROOT}/combined_transfer_rows.tsv" || true
    log "transfer_rows=${OVERNIGHT_ROOT}/combined_transfer_rows.md"
  else
    log "transfer_rows=none"
  fi

  if python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --markdown --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${OVERNIGHT_ROOT}/combined_present_rows.md"; then
    python3 "${HERE}/profile_copy_rows.py" --per-case \
      --top "${PRESENT_ROW_TOP}" --include-zero --sort-by calls \
      --op-prefix ACC_PRESENT \
      "${SUITE_ROOTS[@]}" \
      > "${OVERNIGHT_ROOT}/combined_present_rows.tsv" || true
    log "present_rows=${OVERNIGHT_ROOT}/combined_present_rows.md"
  else
    log "present_rows=none"
  fi
fi
