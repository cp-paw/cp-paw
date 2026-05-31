#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
TEST=${TEST:-si64_bands}
ROOT=$(cd "${HERE}/../../.." && pwd)
NVHPC_STANDARD_ROOT=${NVHPC_STANDARD_ROOT:-"${HERE}/runs/${TEST}-nvhpc-standard-$(date +%Y%m%d-%H%M%S)"}
NSTEPS=${NSTEPS:-3}
EMPTY_BANDS=${EMPTY_BANDS:-1024}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-7200}
GPU_RANKS=${GPU_RANKS:-1}
CPU_RANKS=${CPU_RANKS:-8}
RUN_GPU_ALL=${RUN_GPU_ALL:-no}
GPU_CASES=${GPU_CASES:-"gpu_resident gpu_resident_orthox gpu_resident_addpro_host gpu_resident_pro_host gpu_resident_invbatch_off gpu_resident_no_cusolver"}
CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}
AUTO_BUILD_TARGETS=${AUTO_BUILD_TARGETS:-no}
AUTO_BUILD_JOBS=${AUTO_BUILD_JOBS:-16}

dedup_space_list() {
  local item
  local norm=""
  for item in "$@"; do
    case " ${norm} " in
      *" ${item} "*) ;;
      *) norm="${norm:+${norm} }${item}" ;;
    esac
  done
  echo "${norm}"
}

if [[ ${RUN_GPU_ALL} == yes || ${RUN_GPU_ALL} == true || ${RUN_GPU_ALL} == 1 ]]; then
  GPU_CASES="${GPU_CASES} gpu_all gpu_all_off"
fi

GPU_CASES="${GPU_CASES} gpu_off"
GPU_CASES=$(dedup_space_list ${GPU_CASES})
CPU_CASES=$(dedup_space_list ${CPU_CASES})

REQUIRED_TARGETS=

add_target() {
  local target=$1
  case " ${REQUIRED_TARGETS} " in
    *" ${target} "*) return 0 ;;
  esac
  REQUIRED_TARGETS="${REQUIRED_TARGETS} ${target}"
}

case_target() {
  local case_name=$1
  local ranks=$2
  local target=""
  local gpu_suffix=""

  if (( ranks > 1 )); then
    gpu_suffix="_parallel"
  fi

  case "${case_name}" in
    cpu)
      target="profile${gpu_suffix}"
      ;;
    nvhpc_cpu|nvpl)
      target="nvhpc_profile${gpu_suffix}"
      ;;
    cublas*)
      target="nvhpc_cublas_acc_profile${gpu_suffix}"
      ;;
    cusolver*)
      target="nvhpc_cusolver_acc_profile${gpu_suffix}"
      ;;
    cufftw)
      target="nvhpc_cufftw_profile${gpu_suffix}"
      ;;
    cufft*)
      target="nvhpc_cufft_profile${gpu_suffix}"
      ;;
    nvlamath)
      target="nvhpc_nvlamath_profile${gpu_suffix}"
      ;;
    nvblas)
      target="nvhpc_nvblas_profile${gpu_suffix}"
      ;;
    gpu_resident*)
      target="nvhpc_gpu_acc_residency_profile${gpu_suffix}"
      ;;
    gpu_managed)
      target="nvhpc_gpu_acc_managed_profile${gpu_suffix}"
      ;;
    gpu_unified)
      target="nvhpc_gpu_acc_unified_profile${gpu_suffix}"
      ;;
    gpu_all*)
      target="nvhpc_gpu_all_profile${gpu_suffix}"
      ;;
    gpu*)
      target="nvhpc_gpu_acc_profile${gpu_suffix}"
      ;;
    *)
      echo "unknown" >&2
      return 1
      ;;
  esac

  echo "${target}"
}

target_binary() {
  local target=$1
  case "${target}" in
    profile)
      echo "${ROOT}/bin/profile/paw_profile.x" ;;
    profile_parallel)
      echo "${ROOT}/bin/profile_parallel/ppaw_profile.x" ;;
    nvhpc_profile)
      echo "${ROOT}/bin/nvhpc_profile/paw_nvhpc_profile.x" ;;
    nvhpc_profile_parallel)
      echo "${ROOT}/bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" ;;
    nvhpc_cublas_acc|nvhpc_cublas_acc_profile)
      echo "${ROOT}/bin/nvhpc_cublas_acc_profile/paw_nvhpc_cublas_acc_profile.x" ;;
    nvhpc_cublas_acc_profile_parallel)
      echo "${ROOT}/bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" ;;
    nvhpc_cusolver_acc|nvhpc_cusolver_acc_profile)
      echo "${ROOT}/bin/nvhpc_cusolver_acc_profile/paw_nvhpc_cusolver_acc_profile.x" ;;
    nvhpc_cusolver_acc_profile_parallel)
      echo "${ROOT}/bin/nvhpc_cusolver_acc_profile_parallel/ppaw_nvhpc_cusolver_acc_profile.x" ;;
    nvhpc_cufftw_profile)
      echo "${ROOT}/bin/nvhpc_cufftw_profile/paw_nvhpc_cufftw_profile.x" ;;
    nvhpc_cufftw_profile_parallel)
      echo "${ROOT}/bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" ;;
    nvhpc_cufft_profile)
      echo "${ROOT}/bin/nvhpc_cufft_profile/paw_nvhpc_cufft_profile.x" ;;
    nvhpc_cufft_profile_parallel)
      echo "${ROOT}/bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" ;;
    nvhpc_gpu_acc_profile)
      echo "${ROOT}/bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x" ;;
    nvhpc_gpu_acc_profile_parallel)
      echo "${ROOT}/bin/nvhpc_gpu_acc_profile_parallel/ppaw_nvhpc_gpu_acc_profile.x" ;;
    nvhpc_gpu_acc_residency_profile)
      echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" ;;
    nvhpc_gpu_acc_residency_profile_parallel)
      echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile_parallel/ppaw_nvhpc_gpu_acc_residency_profile.x" ;;
    nvhpc_gpu_all_profile)
      echo "${ROOT}/bin/nvhpc_gpu_all_profile/paw_nvhpc_gpu_all_profile.x" ;;
    nvhpc_gpu_all_profile_parallel)
      echo "${ROOT}/bin/nvhpc_gpu_all_profile_parallel/ppaw_nvhpc_gpu_all_profile.x" ;;
    nvhpc_gpu_acc_managed_profile)
      echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile/paw_nvhpc_gpu_acc_managed_profile.x" ;;
    nvhpc_gpu_acc_managed_profile_parallel)
      echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile_parallel/ppaw_nvhpc_gpu_acc_managed_profile.x" ;;
    nvhpc_gpu_acc_unified_profile)
      echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile/paw_nvhpc_gpu_acc_unified_profile.x" ;;
    nvhpc_gpu_acc_unified_profile_parallel)
      echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile_parallel/ppaw_nvhpc_gpu_acc_unified_profile.x" ;;
    nvhpc_nvlamath_profile)
      echo "${ROOT}/bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" ;;
    nvhpc_nvlamath_profile_parallel)
      echo "${ROOT}/bin/nvhpc_nvlamath_profile_parallel/ppaw_nvhpc_nvlamath_profile.x" ;;
    nvhpc_nvblas_profile)
      echo "${ROOT}/bin/nvhpc_nvblas_profile/paw_nvhpc_nvblas_profile.x" ;;
    nvhpc_nvblas_profile_parallel)
      echo "${ROOT}/bin/nvhpc_nvblas_profile_parallel/ppaw_nvhpc_nvblas_profile.x" ;;
    *)
      echo ""; return 1 ;;
  esac
}

collect_targets() {
  local ranks=$1
  local cases=$2
  local case_name target

  [[ -n ${cases// } ]] || return 0

  for case_name in ${cases}; do
    if ! target=$(case_target "${case_name}" "${ranks}"); then
      echo "Unknown case in GPU/CPU selection: ${case_name}" >&2
      exit 1
    fi
    add_target "${target}"
  done
}

ensure_binaries() {
  local build_targets
  local target exe missing=()

  if [[ -z "${REQUIRED_TARGETS}" ]]; then
    echo "No benchmark targets were discovered from selected cases." >&2
    exit 1
  fi

  for target in ${REQUIRED_TARGETS}; do
    exe=$(target_binary "${target}") || exit 1
    if [[ ! -x "${exe}" ]]; then
      missing+=("${target}")
    fi
  done

  if (( ${#missing[@]} == 0 )); then
    return
  fi

  if [[ ${AUTO_BUILD_TARGETS} == yes || ${AUTO_BUILD_TARGETS} == true || ${AUTO_BUILD_TARGETS} == 1 ]]; then
    build_targets=$(printf '%s\n' "${missing[@]}" | sort -u)
    log "Auto-building required profile targets: ${build_targets}"
    for target in ${build_targets}; do
      log "Building target=${target}"
      (cd "${ROOT}" && CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -z -j "${AUTO_BUILD_JOBS}" -c "${target}")
    done
    return
  fi

  log "Missing required binaries for selected cases:"
  for target in "${missing[@]}"; do
    log "  target=${target} exe=$(target_binary "${target}")"
  done
  log "Set AUTO_BUILD_TARGETS=yes to build them, then rerun."
  exit 1
}

mkdir -p "${NVHPC_STANDARD_ROOT}"
echo "${NVHPC_STANDARD_ROOT}" > "${HERE}/runs/latest_nvhpc_standard"

COMBINED="${NVHPC_STANDARD_ROOT}/combined_benchmark.tsv"
LOG="${NVHPC_STANDARD_ROOT}/nvhpc_standard.log"
: > "${LOG}"

log() {
  printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${LOG}"
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
  local root="${NVHPC_STANDARD_ROOT}/${suite}"
  local suite_log="${NVHPC_STANDARD_ROOT}/${suite}.log"

  if [[ -z ${cases// } ]]; then
    log "SKIP  suite=${suite} empty case list"
    echo "skipped" > "${root}.status"
    return 0
  fi

  log "START suite=${suite} empty_bands=${EMPTY_BANDS} nsteps=${NSTEPS} ranks=${ranks} cases=${cases}"
  if env TEST="${TEST}" EMPTY_BANDS="${EMPTY_BANDS}" NSTEPS="${NSTEPS}" \
      RANKS="${ranks}" REPEATS="${REPEATS}" CASES="${cases}" TIMEOUT="${TIMEOUT}" \
      REQUIRE_CASES="yes" RUN_ROOT="${root}" "${HERE}/run_benchmark.sh" > "${suite_log}" 2>&1; then
    log "DONE  suite=${suite}"
    echo 0 > "${root}.status"
  else
    local status=$?
    log "FAIL  suite=${suite} status=${status}"
    echo "${status}" > "${root}.status"
  fi
  append_suite "${suite}" "${root}/benchmark.tsv"
}

collect_targets "${GPU_RANKS}" "${GPU_CASES}"
collect_targets "1" "${CPU_CASES}"
collect_targets "${CPU_RANKS}" "${CPU_CASES}"
ensure_binaries

run_suite "gpu_${GPU_RANKS}rank" "${GPU_RANKS}" "${GPU_CASES}"
run_suite "cpu_1rank" 1 "${CPU_CASES}"
run_suite "cpu_${CPU_RANKS}rank_ref" "${CPU_RANKS}" "${CPU_CASES}"

log "ALL DONE root=${NVHPC_STANDARD_ROOT}"
if [[ -f "${COMBINED}" ]]; then
  python3 "${HERE}/benchmark_markdown.py" "${COMBINED}" \
    > "${NVHPC_STANDARD_ROOT}/combined_benchmark.md" || true
  log "combined=${COMBINED}"
fi
