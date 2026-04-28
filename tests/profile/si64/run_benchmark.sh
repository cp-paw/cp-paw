#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
TEST=${TEST:-si64}
NSTEPS=${NSTEPS:-20}
RANKS=${RANKS:-4}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-1800}
RUN_ROOT=${RUN_ROOT:-"${HERE}/runs/${TEST}-nstep${NSTEPS}-${RANKS}ranks-$(date +%Y%m%d-%H%M%S)"}
MPI_ARGS=${MPI_ARGS:---mca coll ^hcoll}
CASES=${CASES:-"nvpl cublas cublas_conservative cublas_off"}

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

MPIRUN=${MPIRUN:-$(default_mpirun)}

serial_exe() {
  case "$1" in
    nvpl) echo "${ROOT}/bin/nvhpc_profile/paw_nvhpc_profile.x" ;;
    nvblas) echo "${ROOT}/bin/nvhpc_nvblas_profile/paw_nvhpc_nvblas_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile/paw_nvhpc_cufftw_profile.x" ;;
    cufft|cufft_off) echo "${ROOT}/bin/nvhpc_cufft_profile/paw_nvhpc_cufft_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_cufft_cublas_acc_profile/paw_nvhpc_cufft_cublas_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile/paw_nvhpc_cublas_acc_profile.x" ;;
    *) echo "unknown case $1" >&2; return 1 ;;
  esac
}

parallel_exe() {
  case "$1" in
    nvpl) echo "${ROOT}/bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" ;;
    nvblas) echo "${ROOT}/bin/nvhpc_nvblas_profile_parallel/ppaw_nvhpc_nvblas_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" ;;
    cufft|cufft_off) echo "${ROOT}/bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_cufft_cublas_acc_profile_parallel/ppaw_nvhpc_cufft_cublas_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" ;;
    *) echo "unknown case $1" >&2; return 1 ;;
  esac
}

run_command() {
  local case_name=$1
  local exe=$2
  local cmd
  if [[ "${RANKS}" -gt 1 ]]; then
    # shellcheck disable=SC2086
    cmd="${MPIRUN} ${MPI_ARGS} -np ${RANKS} ${exe} ${TEST}.cntl"
  else
    cmd="${exe} ${TEST}.cntl"
  fi
  if [[ "${case_name}" = nvblas ]]; then
    echo "$(dirname "${exe}")/paw_nvblas.sh ${cmd}"
  else
    echo "${cmd}"
  fi
}

case_env() {
  case "$1" in
    cublas) echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_ACC_MINFLOP:-1e7}" ;;
    cublas_conservative) echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}" ;;
    cublas_off) echo "CPPAW_CUBLAS_ACC=0" ;;
    cufft) echo "CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0}" ;;
    cufft_off) echo "CPPAW_CUFFT_ACC=0" ;;
    gpu) echo "CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0} CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_ACC_MINFLOP:-1e7}" ;;
    gpu_conservative) echo "CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0} CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}" ;;
    gpu_no_cufft) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_ACC_MINFLOP:-1e7}" ;;
    gpu_no_cublas) echo "CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0} CPPAW_CUBLAS_ACC=0" ;;
    gpu_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0" ;;
    *) echo "" ;;
  esac
}

prepare_case() {
  local dir=$1
  mkdir -p "${dir}"
  cp "${HERE}/${TEST}.cntl" "${HERE}/${TEST}.strc" \
     "${HERE}/profile_summary.py" "${HERE}/benchmark_summary.py" "${dir}/"
  cp "${ROOT}/tests/fulltests/si2/stp.cntl" "${dir}/"
  perl -0pi -e "s/NSTEP\\s*=\\s*\\d+/NSTEP=${NSTEPS}/" "${dir}/${TEST}.cntl"
}

mkdir -p "${RUN_ROOT}"
echo "${RUN_ROOT}" > "${HERE}/runs/latest"

for case_name in ${CASES}; do
  exe=$(if [[ "${RANKS}" -gt 1 ]]; then parallel_exe "${case_name}"; else serial_exe "${case_name}"; fi)
  if [[ ! -x "${exe}" ]]; then
    echo "Skipping ${case_name}: executable not found: ${exe}" >&2
    continue
  fi
  for repeat in $(seq 1 "${REPEATS}"); do
    run_dir="${RUN_ROOT}/${case_name}/rep$(printf "%02d" "${repeat}")"
    prepare_case "${run_dir}"
    (
      cd "${run_dir}"
      env_line=$(case_env "${case_name}")
      {
        echo "case=${case_name}"
        echo "repeat=${repeat}"
        echo "nsteps=${NSTEPS}"
        echo "ranks=${RANKS}"
        echo "exe=${exe}"
        echo "env=${env_line}"
        echo "start=$(date -Is)"
      } > run.env
      cmd=$(run_command "${case_name}" "${exe}")
      echo "running ${case_name} repeat ${repeat}: ${cmd}"
      if [[ -n "${env_line}" ]]; then
        # shellcheck disable=SC2086
        /usr/bin/time -p env CPPAW_ACCEL_PROFILE_FILE="${case_name}_profile" ${env_line} \
          timeout "${TIMEOUT}s" ${cmd} > out.log 2> err.log
      else
        # shellcheck disable=SC2086
        /usr/bin/time -p env CPPAW_ACCEL_PROFILE_FILE="${case_name}_profile" \
          timeout "${TIMEOUT}s" ${cmd} > out.log 2> err.log
      fi
      python3 profile_summary.py "${case_name}_profile"*.csv > summary.txt
      grep -E "^real |^user |^sys " err.log > time.txt || true
      tail -40 out.log > out.tail.txt || true
      tail -80 err.log > err.tail.txt || true
      echo "end=$(date -Is)" >> run.env
    )
  done
done

python3 "${HERE}/benchmark_summary.py" "${RUN_ROOT}" | tee "${RUN_ROOT}/benchmark.tsv"
echo "Benchmark data: ${RUN_ROOT}"
