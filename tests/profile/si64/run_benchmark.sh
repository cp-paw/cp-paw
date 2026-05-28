#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)
TEST=${TEST:-si64}
CNTL_FILE=${CNTL_FILE:-"${HERE}/${TEST}.cntl"}
STRC_FILE=${STRC_FILE:-"${HERE}/${TEST}.strc"}
if [[ ! -f "${STRC_FILE}" && -f "${HERE}/si64.strc" ]]; then
  STRC_FILE="${HERE}/si64.strc"
fi
NSTEPS=${NSTEPS:-20}
RANKS=${RANKS:-1}
REPEATS=${REPEATS:-1}
TIMEOUT=${TIMEOUT:-1800}
RUN_ROOT=${RUN_ROOT:-"${HERE}/runs/${TEST}-nstep${NSTEPS}-${RANKS}ranks-$(date +%Y%m%d-%H%M%S)"}
MPI_ARGS=${MPI_ARGS:---mca coll ^hcoll}
CASES=${CASES:-"cpu nvhpc_cpu gpu_resident gpu_off"}
TIMEOUT_PREFIX=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_PREFIX="timeout ${TIMEOUT}s"
fi

TIME_CMD=${TIME_CMD:-}
if [[ -z "${TIME_CMD}" ]]; then
  if [[ -x /usr/bin/time ]]; then
    TIME_CMD=/usr/bin/time
  else
    TIME_CMD=$(type -P time || true)
  fi
fi
if [[ -z "${TIME_CMD}" ]]; then
  echo "External time command not found; set TIME_CMD or install GNU time." >&2
  exit 1
fi

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

system_mpirun() {
  command -v mpirun 2>/dev/null || echo mpirun
}

CPU_MPIRUN=${CPU_MPIRUN:-$(system_mpirun)}
MPIRUN=${MPIRUN:-$(default_mpirun)}

case_mpirun() {
  case "$1" in
    cpu) echo "${CPU_MPIRUN}" ;;
    *) echo "${MPIRUN}" ;;
  esac
}

serial_exe() {
  case "$1" in
    cpu) echo "${ROOT}/bin/profile/paw_profile.x" ;;
    nvhpc_cpu|nvpl) echo "${ROOT}/bin/nvhpc_profile/paw_nvhpc_profile.x" ;;
    nvblas) echo "${ROOT}/bin/nvhpc_nvblas_profile/paw_nvhpc_nvblas_profile.x" ;;
    nvlamath) echo "${ROOT}/bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile/paw_nvhpc_cufftw_profile.x" ;;
    cufft|cufft_off) echo "${ROOT}/bin/nvhpc_cufft_profile/paw_nvhpc_cufft_profile.x" ;;
    gpu_all*) echo "${ROOT}/bin/nvhpc_gpu_all_profile/paw_nvhpc_gpu_all_profile.x" ;;
    gpu_resident*) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" ;;
    gpu_managed*) echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile/paw_nvhpc_gpu_acc_managed_profile.x" ;;
    gpu_unified*) echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile/paw_nvhpc_gpu_acc_unified_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile/paw_nvhpc_cublas_acc_profile.x" ;;
    cusolver*) echo "${ROOT}/bin/nvhpc_cusolver_acc_profile/paw_nvhpc_cusolver_acc_profile.x" ;;
    *) echo "unknown case $1" >&2; return 1 ;;
  esac
}

parallel_exe() {
  case "$1" in
    cpu) echo "${ROOT}/bin/profile_parallel/ppaw_profile.x" ;;
    nvhpc_cpu|nvpl) echo "${ROOT}/bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" ;;
    nvblas) echo "${ROOT}/bin/nvhpc_nvblas_profile_parallel/ppaw_nvhpc_nvblas_profile.x" ;;
    nvlamath) echo "${ROOT}/bin/nvhpc_nvlamath_profile_parallel/ppaw_nvhpc_nvlamath_profile.x" ;;
    cufftw) echo "${ROOT}/bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" ;;
    cufft|cufft_off) echo "${ROOT}/bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" ;;
    gpu_all*) echo "${ROOT}/bin/nvhpc_gpu_all_profile_parallel/ppaw_nvhpc_gpu_all_profile.x" ;;
    gpu_resident*) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile_parallel/ppaw_nvhpc_gpu_acc_residency_profile.x" ;;
    gpu_managed*) echo "${ROOT}/bin/nvhpc_gpu_acc_managed_profile_parallel/ppaw_nvhpc_gpu_acc_managed_profile.x" ;;
    gpu_unified*) echo "${ROOT}/bin/nvhpc_gpu_acc_unified_profile_parallel/ppaw_nvhpc_gpu_acc_unified_profile.x" ;;
    gpu*) echo "${ROOT}/bin/nvhpc_gpu_acc_profile_parallel/ppaw_nvhpc_gpu_acc_profile.x" ;;
    cublas*) echo "${ROOT}/bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" ;;
    cusolver*) echo "${ROOT}/bin/nvhpc_cusolver_acc_profile_parallel/ppaw_nvhpc_cusolver_acc_profile.x" ;;
    *) echo "unknown case $1" >&2; return 1 ;;
  esac
}

case_note() {
  case "$1" in
    nvhpc_cpu)
      echo "NVIDIA HPC SDK CPU build; actual CPU BLAS/LAPACK/FFT libraries are listed below."
      ;;
    nvpl)
      echo "Legacy alias for nvhpc_cpu; on x86 this may use OpenBLAS/FFTW fallback rather than NVPL."
      ;;
    gpu_resident)
      echo "Recommended one-GPU NVHPC profile path with OpenACC residency enabled."
      ;;
    gpu_resident_nosync)
      echo "Residency diagnostic that disables the explicit post-cuBLAS device synchronization."
      ;;
    *_projection_conservative)
      echo "cuBLAS diagnostic: raises only the projection GEMM offload threshold."
      ;;
    *_overlap_conservative)
      echo "cuBLAS diagnostic: raises only the overlap/orthogonalization offload threshold."
      ;;
    *_addproduct_conservative)
      echo "cuBLAS diagnostic: raises only the additive product offload threshold."
      ;;
    *_matmul_conservative)
      echo "cuBLAS diagnostic: raises only the generic MATMUL offload threshold."
      ;;
    gpu_all*)
      echo "All-library diagnostic build; includes cuFFTW/NVLAMATH and is not the recommended default."
      ;;
    *)
      echo ""
      ;;
  esac
}

run_command() {
  local case_name=$1
  local exe=$2
  local cmd mpi_run
  if [[ "${RANKS}" -gt 1 ]]; then
    mpi_run=$(case_mpirun "${case_name}")
    # shellcheck disable=SC2086
    cmd="${mpi_run} ${MPI_ARGS} -np ${RANKS} ${exe} ${TEST}.cntl"
  else
    cmd="${exe} ${TEST}.cntl"
  fi
  if [[ "${case_name}" = nvblas ]]; then
    echo "$(dirname "${exe}")/paw_nvblas.sh ${cmd}"
  else
    echo "${cmd}"
  fi
}

cusolver_env() {
  local min_n=$1
  local env_line="CPPAW_CUSOLVER_ACC_MIN_N=${min_n}"
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK=${CPPAW_CUSOLVER_ACC_CHECK}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK_TOL:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK_TOL=${CPPAW_CUSOLVER_ACC_CHECK_TOL}"
  fi
  echo "${env_line}"
}

cufft_min_default() {
  echo "${CPPAW_CUFFT_CONSERVATIVE_MIN_ELEMENTS:-1000000}"
}

cufft_env() {
  local min_elements
  min_elements=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-$(cufft_min_default)}
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${min_elements}"
}

cufft_force_env() {
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-0}"
}

cufft3d_env() {
  local min_elements min_3d_elements
  min_elements=${CPPAW_CUFFT_ACC_MIN_ELEMENTS:-$(cufft_min_default)}
  min_3d_elements=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-$(cufft_min_default)}
  echo "CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_MIN_ELEMENTS=${min_elements} CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${min_3d_elements}"
}

cublas_env() {
  echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_ACC_MINFLOP:-1e7}"
}

cublas_conservative_env() {
  echo "CPPAW_CUBLAS_ACC_MINFLOP=${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}"
}

cublas_projection_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_PROJECTION_MINFLOP=${CPPAW_CUBLAS_PROJECTION_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_overlap_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_OVERLAP_MINFLOP=${CPPAW_CUBLAS_OVERLAP_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_addproduct_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP=${CPPAW_CUBLAS_ADDPRODUCT_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

cublas_matmul_conservative_env() {
  echo "$(cublas_env) CPPAW_CUBLAS_ACC_MATMUL_MINFLOP=${CPPAW_CUBLAS_MATMUL_CONSERVATIVE_MINFLOP:-${CPPAW_CUBLAS_CONSERVATIVE_MINFLOP:-1e8}}"
}

case_env() {
  case "$1" in
    cublas) cublas_env ;;
    cublas_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    cublas_conservative) cublas_conservative_env ;;
    cublas_projection_conservative) cublas_projection_conservative_env ;;
    cublas_overlap_conservative) cublas_overlap_conservative_env ;;
    cublas_addproduct_conservative) cublas_addproduct_conservative_env ;;
    cublas_matmul_conservative) cublas_matmul_conservative_env ;;
    cublas_off) echo "CPPAW_CUBLAS_ACC=0" ;;
    cusolver) cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}" ;;
    cusolver_conservative) cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}" ;;
    cusolver_off) echo "CPPAW_CUSOLVER_ACC=0" ;;
    cufft) cufft_env ;;
    cufft_force_all) cufft_force_env ;;
    cufft_off) echo "CPPAW_CUFFT_ACC=0" ;;
    gpu) cublas_env ;;
    gpu_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_projection_conservative) cublas_projection_conservative_env ;;
    gpu_overlap_conservative) cublas_overlap_conservative_env ;;
    gpu_addproduct_conservative) cublas_addproduct_conservative_env ;;
    gpu_matmul_conservative) cublas_matmul_conservative_env ;;
    gpu_force_all) echo "$(cufft_force_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_conservative) echo "$(cublas_conservative_env) $(cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}")" ;;
    gpu_all) echo "$(cufft_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_nosync) echo "$(cufft_env) $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_resident) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_nosync) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_resident_projection_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_projection_conservative_env)" ;;
    gpu_resident_overlap_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_overlap_conservative_env)" ;;
    gpu_resident_addproduct_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_addproduct_conservative_env)" ;;
    gpu_resident_matmul_conservative) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_matmul_conservative_env)" ;;
    gpu_resident_force_all) echo "CPPAW_GPU_RESIDENCY=1 $(cufft_force_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_resident_off) echo "CPPAW_GPU_RESIDENCY=0 CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_managed|gpu_unified) cublas_env ;;
    gpu_no_cufft) echo "CPPAW_CUFFT_ACC=0 $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_no_cublas) echo "$(cufft_force_env) CPPAW_CUBLAS_ACC=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_no_cusolver) echo "$(cufft_force_env) $(cublas_env) CPPAW_CUSOLVER_ACC=0" ;;
    gpu_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    *) echo "" ;;
  esac
}

iso_now() {
  date -u +%Y-%m-%dT%H:%M:%SZ
}

prepare_case() {
  local dir=$1
  mkdir -p "${dir}"
  if [[ ! -f "${CNTL_FILE}" ]]; then
    echo "Control file not found: ${CNTL_FILE}" >&2
    return 1
  fi
  if [[ ! -f "${STRC_FILE}" ]]; then
    echo "Structure file not found: ${STRC_FILE}" >&2
    return 1
  fi
  cp "${CNTL_FILE}" "${dir}/${TEST}.cntl"
  cp "${STRC_FILE}" "${dir}/${TEST}.strc"
  cp "${HERE}/profile_summary.py" "${HERE}/benchmark_summary.py" "${dir}/"
  cp "${ROOT}/tests/fulltests/si2/stp.cntl" "${dir}/"
  perl -0pi -e "s/NSTEP\\s*=\\s*\\d+/NSTEP=${NSTEPS}/" "${dir}/${TEST}.cntl"
  if [[ -n "${EMPTY_BANDS:-}" ]]; then
    perl -0pi -e "s/EMPTY\\s*=\\s*\\d+/EMPTY=${EMPTY_BANDS}/" "${dir}/${TEST}.strc"
  fi
}

capture_metadata() {
  {
    echo "date=$(iso_now)"
    echo "hostname=$(hostname)"
    echo "test=${TEST}"
    echo "nsteps=${NSTEPS}"
    echo "ranks=${RANKS}"
    echo "cases=${CASES}"
    echo "threads=OMP_NUM_THREADS=${OMP_NUM_THREADS} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS} BLIS_NUM_THREADS=${BLIS_NUM_THREADS} VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS} NVPL_NUM_THREADS=${NVPL_NUM_THREADS}"
    echo
    uname -a
    echo
    for case_name in ${CASES}; do
      exe=$(if [[ "${RANKS}" -gt 1 ]]; then parallel_exe "${case_name}"; else serial_exe "${case_name}"; fi)
      echo "case=${case_name}"
      echo "exe=${exe}"
      note=$(case_note "${case_name}")
      [[ -n "${note}" ]] && echo "note=${note}"
      if [[ -x "${exe}" ]]; then
        ldd "${exe}" 2>/dev/null | grep -E "blas|lapack|fftw|cufft|cusolver|cublas|nvpl|openblas" || true
      else
        echo "missing"
      fi
      echo
    done
  } > "${RUN_ROOT}/metadata.txt" 2>&1
}

mkdir -p "${RUN_ROOT}"
echo "${RUN_ROOT}" > "${HERE}/runs/latest"
capture_metadata

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
        echo "threads=OMP_NUM_THREADS=${OMP_NUM_THREADS} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS} BLIS_NUM_THREADS=${BLIS_NUM_THREADS} VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS} NVPL_NUM_THREADS=${NVPL_NUM_THREADS}"
        echo "start=$(iso_now)"
      } > run.env
      cmd=$(run_command "${case_name}" "${exe}")
      echo "running ${case_name} repeat ${repeat}: ${cmd}"
      if [[ -n "${env_line}" ]]; then
        # shellcheck disable=SC2086
        "${TIME_CMD}" -p env CPPAW_ACCEL_PROFILE_FILE="${case_name}_profile" ${env_line} \
          ${TIMEOUT_PREFIX} ${cmd} > out.log 2> err.log
      else
        # shellcheck disable=SC2086
        "${TIME_CMD}" -p env CPPAW_ACCEL_PROFILE_FILE="${case_name}_profile" \
          ${TIMEOUT_PREFIX} ${cmd} > out.log 2> err.log
      fi
      python3 profile_summary.py "${case_name}_profile"*.csv > summary.txt
      tr -d '\000' < err.log | grep -E "^real |^user |^sys " > time.txt || true
      tail -40 out.log > out.tail.txt || true
      tail -80 err.log > err.tail.txt || true
      echo "end=$(iso_now)" >> run.env
    )
  done
done

python3 "${HERE}/benchmark_summary.py" "${RUN_ROOT}" | tee "${RUN_ROOT}/benchmark.tsv"
python3 "${HERE}/benchmark_markdown.py" "${RUN_ROOT}/benchmark.tsv" \
  > "${RUN_ROOT}/benchmark.md" || true
echo "Benchmark data: ${RUN_ROOT}"
