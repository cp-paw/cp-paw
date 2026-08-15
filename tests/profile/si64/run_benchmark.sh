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
REQUIRE_CASES=${REQUIRE_CASES:-no}
DRY_RUN=${DRY_RUN:-no}
RUN_ROOT=${RUN_ROOT:-"${HERE}/runs/${TEST}-nstep${NSTEPS}-${RANKS}ranks-$(date +%Y%m%d-%H%M%S)"}
MPI_ARGS=${MPI_ARGS:---mca coll ^hcoll}
CASES=${CASES:-"cpu nvhpc_cpu gpu_resident gpu_off"}
EXPECTED_ENERGY_USER_SET=${EXPECTED_ENERGY+x}
EXPECTED_ENERGY=${EXPECTED_ENERGY:-}
ENERGY_TOL=${ENERGY_TOL:-1e-5}
ENERGY_CHECK=${ENERGY_CHECK:-yes}
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
if [[ -z "${TIME_CMD}" ]] && command -v python3 >/dev/null 2>&1; then
  mkdir -p "${RUN_ROOT}"
  TIME_CMD="$(cd "${RUN_ROOT}" && pwd)/time_posix.py"
  cat > "${TIME_CMD}" <<'PY'
#!/usr/bin/env python3
import resource
import subprocess
import sys
import time

args = sys.argv[1:]
if args and args[0] == "-p":
    args = args[1:]

start_time = time.time()
start_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
proc = subprocess.run(args)
end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
elapsed = time.time() - start_time
user = end_usage.ru_utime - start_usage.ru_utime
system = end_usage.ru_stime - start_usage.ru_stime

sys.stderr.write(f"real {elapsed:.2f}\n")
sys.stderr.write(f"user {user:.2f}\n")
sys.stderr.write(f"sys {system:.2f}\n")
raise SystemExit(proc.returncode)
PY
  chmod +x "${TIME_CMD}"
fi
if [[ -z "${TIME_CMD}" ]]; then
  case "${DRY_RUN}" in
    yes|true|1) ;;
    *)
      echo "External time command not found; set TIME_CMD or install GNU time." >&2
      exit 1
      ;;
  esac
fi

case "${EXPECTED_ENERGY}" in
  none|off|skip)
    EXPECTED_ENERGY=
    ENERGY_CHECK=no
    ;;
esac

case "${ENERGY_CHECK}" in
  no|false|0) EXPECTED_ENERGY= ;;
  *)
    case "${TEST}" in
      si64|si64_bands)
        if [[ -z "${EXPECTED_ENERGY_USER_SET}" && "${NSTEPS}" != "1" ]]; then
          EXPECTED_ENERGY=
        else
          EXPECTED_ENERGY=${EXPECTED_ENERGY:-302.280854}
        fi
        ;;
    esac
    ;;
esac

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

system_mpirun() {
  local candidate
  for candidate in \
      "${CONDA_PREFIX:-}/bin/mpirun" \
      "${MAMBA_ROOT_PREFIX:-}/envs/cppaw-gccmpi/bin/mpirun" \
      "${HOME:-}/micromamba/envs/cppaw-gccmpi/bin/mpirun" \
      /usr/bin/mpirun \
      /usr/local/bin/mpirun; do
    if [[ -x "${candidate}" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  command -v mpirun 2>/dev/null || echo mpirun
}

CPU_MPIRUN_USER_SET=${CPU_MPIRUN+x}
CPU_MPI_LIBDIR_USER_SET=${CPU_MPI_LIBDIR+x}
CPU_MPIRUN=${CPU_MPIRUN:-$(system_mpirun)}
MPIRUN=${MPIRUN:-$(default_mpirun)}

mpirun_libdir() {
  local mpirun=$1
  local libdir
  if [[ "${mpirun}" = */bin/mpirun ]]; then
    libdir=$(cd "$(dirname "${mpirun}")/../lib" 2>/dev/null && pwd || true)
    if [[ -n "${libdir}" && -d "${libdir}" ]]; then
      echo "${libdir}"
    fi
  fi
}

CPU_MPI_LIBDIR=${CPU_MPI_LIBDIR:-$(mpirun_libdir "${CPU_MPIRUN}")}

exe_mpi_libdir() {
  local exe=$1
  local mpi_lib
  [[ -x "${exe}" ]] || return 0
  mpi_lib=$(ldd "${exe}" 2>/dev/null | awk '/libmpi[.]so/ { print $3; exit }')
  if [[ -n "${mpi_lib}" && -f "${mpi_lib}" ]]; then
    dirname "${mpi_lib}"
  fi
}

mpirun_from_libdir() {
  local libdir=$1
  local fallback=$2
  local prefix candidate
  if [[ -n "${libdir}" ]]; then
    prefix=${libdir%/lib}
    for candidate in \
        "${prefix}/bin/mpirun" \
        "${prefix}/ompi/bin/mpirun"; do
      if [[ -x "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    done
  fi
  echo "${fallback}"
}

cpu_mpi_libdir_for_exe() {
  local exe=$1
  if [[ -n "${CPU_MPI_LIBDIR_USER_SET}" || -n "${CPU_MPIRUN_USER_SET}" ]]; then
    echo "${CPU_MPI_LIBDIR}"
    return 0
  fi
  exe_mpi_libdir "${exe}"
}

cpu_mpirun_for_exe() {
  local exe=$1
  local libdir
  if [[ -n "${CPU_MPIRUN_USER_SET}" ]]; then
    echo "${CPU_MPIRUN}"
    return 0
  fi
  libdir=$(cpu_mpi_libdir_for_exe "${exe}")
  mpirun_from_libdir "${libdir}" "${CPU_MPIRUN}"
}

case_mpirun() {
  case "$1" in
    cpu) cpu_mpirun_for_exe "${2:-}" ;;
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
    gpu_resident*|gpu_psim_propagate|gpu_hpsi_psim_propagate) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" ;;
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
    gpu_resident*|gpu_psim_propagate|gpu_hpsi_psim_propagate) echo "${ROOT}/bin/nvhpc_gpu_acc_residency_profile_parallel/ppaw_nvhpc_gpu_acc_residency_profile.x" ;;
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
    gpu_resident_1coverlap)
      echo "Residency diagnostic that explicitly enables one-center overlap cuBLAS offload."
      ;;
    gpu_resident_1coverlap_batch)
      echo "Opt-in diagnostic that batches the three orthogonalization one-center overlap matrices in one cuBLAS path."
      ;;
    gpu_resident_addoproj)
      echo "Opt-in diagnostic that offloads orthogonalization WAVES_ADDOPROJ projector updates through cuBLAS slice GEMMs."
      ;;
    gpu_resident_nosync)
      echo "Residency diagnostic that disables the explicit post-cuBLAS device synchronization."
      ;;
    gpu_resident_pro_host)
      echo "Residency diagnostic that keeps projector expansion on host while retaining resident cuBLAS projection."
      ;;
    gpu_resident_proj)
      echo "Residency diagnostic that keeps THIS%PROJ present after WAVES\$PROJECTIONS for downstream DENMAT/ADDPRO/off-site reuse."
      ;;
    gpu_resident_addpro_host)
      echo "Residency diagnostic that keeps the GPU projector cache for projections but disables its WAVES_ADDPRO reuse."
      ;;
    gpu_resident_addpro_hpsi_host)
      echo "Residency diagnostic that disables only HPSI-side WAVES_ADDPRO cache reuse."
      ;;
    gpu_resident_addpro_opsi_host)
      echo "Residency diagnostic that disables only OPSI-side WAVES_ADDPRO cache reuse."
      ;;
    gpu_resident_opsi_addpro_host)
      echo "OPSI-residency diagnostic with OPSI-side WAVES_ADDPRO cache reuse disabled."
      ;;
    gpu_resident_forcepsi_host)
      echo "Residency diagnostic that disables the force-loop PSI0 resident data region."
      ;;
    gpu_resident_orthoconst)
      echo "Residency diagnostic that enables resident CHICHI/U inputs in WAVES_ORTHO_X."
      ;;
    gpu_resident_orthox)
      echo "Residency default that keeps the WAVES_ORTHO_X iteration workspace on the GPU."
      ;;
    gpu_resident_orthox_off)
      echo "Residency diagnostic that disables the WAVES_ORTHO_X iteration workspace residency."
      ;;
    gpu_resident_orthox_nosync)
      echo "Residency diagnostic that combines WAVES_ORTHO_X workspace residency with disabled post-cuBLAS device synchronization."
      ;;
    gpu_resident_opsi)
      echo "Residency diagnostic that keeps the orthogonalization OPSI wavefunction on the GPU from its build through ADDOPSI."
      ;;
    gpu_resident_hpsi)
      echo "Residency diagnostic that keeps HPSI on the GPU from WAVES_ADDPRO through the immediate expectation/Hamiltonian overlaps."
      ;;
    gpu_resident_denmat_energy)
      echo "Opt-in diagnostic that offloads the time-inversion one-center DENMAT energy/Lambda contraction with OpenACC."
      ;;
    gpu_resident_hpsi_denmat_energy)
      echo "Combined HPSI residency plus DENMAT energy/Lambda OpenACC diagnostic."
      ;;
    gpu_resident_denmat_energy_offden_blas)
      echo "Combined DENMAT energy/Lambda OpenACC diagnostic plus scalar TINV off-site DENMAT BLAS prototype."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_blas)
      echo "Combined HPSI residency, DENMAT energy/Lambda OpenACC diagnostic, and scalar TINV off-site DENMAT BLAS prototype."
      ;;
    gpu_resident_offden_cublas)
      echo "Residency diagnostic that tries cuBLAS for the scalar TINV off-site DENMAT BLAS prototype."
      ;;
    gpu_resident_hpsi_offden_cublas)
      echo "Combined HPSI residency plus scalar TINV off-site DENMAT cuBLAS diagnostic."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas)
      echo "Combined HPSI residency, DENMAT energy/Lambda diagnostic, and scalar TINV off-site DENMAT cuBLAS diagnostic."
      ;;
    gpu_resident_hpsi_offden_cublas_batch)
      echo "Combined HPSI residency plus stacked cuBLAS diagnostic for scalar TINV off-site DENMAT."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_batch)
      echo "Combined HPSI residency, DENMAT energy/Lambda diagnostic, and stacked cuBLAS off-site DENMAT diagnostic."
      ;;
    gpu_resident_hpsi_offden_cublas_devicepack)
      echo "Combined HPSI residency plus stacked cuBLAS off-site DENMAT with OpenACC device packing."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack)
      echo "Combined HPSI residency, DENMAT energy/Lambda diagnostic, and stacked cuBLAS off-site DENMAT with OpenACC device packing."
      ;;
    gpu_resident_hpsi_offden_cublas_devicepack_accum)
      echo "Device-pack off-site DENMAT diagnostic that accumulates the stacked cuBLAS result into a real matrix pack on the GPU."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum)
      echo "Combined DENMAT energy/device-pack diagnostic with GPU-side off-site real-matrix accumulation."
      ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj)
      echo "Device-pack off-site DENMAT diagnostic with HPSI and persistent THIS%PROJ residency enabled."
      ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj_accum)
      echo "Device-pack off-site DENMAT diagnostic with HPSI, persistent THIS%PROJ, and GPU-side real-matrix accumulation."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj)
      echo "Combined DENMAT energy/device-pack off-site diagnostic with persistent THIS%PROJ residency enabled."
      ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum)
      echo "Combined DENMAT energy/device-pack off-site diagnostic with persistent THIS%PROJ and GPU-side real-matrix accumulation."
      ;;
    gpu_resident_offden_blas)
      echo "Residency diagnostic that enables the opt-in scalar TINV off-site DENMAT BLAS prototype."
      ;;
    gpu_resident_hpsi_offden_blas)
      echo "Combined HPSI residency plus scalar TINV off-site DENMAT BLAS prototype."
      ;;
    gpu_psim_propagate|gpu_resident_psim)
      echo "Diagnostic that propagates PSIM on the GPU and copies it back before orthogonalization."
      ;;
    gpu_hpsi_psim_propagate|gpu_resident_hpsi_psim)
      echo "Diagnostic that combines HPSI residency with GPU PSIM propagation."
      ;;
    gpu_resident_psim_phase)
      echo "Diagnostic that leaves the propagated PSIM present until the following orthogonalization copies it back."
      ;;
    gpu_resident_hpsi_psim_phase)
      echo "Diagnostic that combines HPSI residency with cross-phase PSIM propagation residency."
      ;;
    gpu_resident_hpsi_opsi_psim_phase)
      echo "Diagnostic that combines HPSI, OPSI, and cross-phase PSIM propagation residency."
      ;;
    gpu_resident_hpsi_opsi)
      echo "Residency diagnostic that combines HPSI and OPSI wavefunction residency switches."
      ;;
    gpu_resident_hpsi_opsi_proj)
      echo "Residency diagnostic that combines HPSI, OPSI, and persistent THIS%PROJ residency."
      ;;
    gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum)
      echo "Combined HPSI/OPSI residency plus off-site DENMAT device packing and GPU-side real-matrix accumulation."
      ;;
    gpu_resident_stack)
      echo "Focused residency stack keyword: setup PSI0, PSI0-to-PRINFO, HPSI/OPSI, PSIM switch, PROJ, DENMAT energy, and off-site device-pack accumulation."
      ;;
    gpu_resident_stack_psi0_ortho_host)
      echo "Focused residency stack with cross-orthogonalization PSI0 residency disabled for copy-boundary diagnostics."
      ;;
    gpu_resident_stack_psi0_prinfo_host)
      echo "Focused residency stack with PSI0-to-PRINFO residency disabled for WRITEPDOS projection diagnostics."
      ;;
    gpu_resident_stack_setup_host)
      echo "Focused residency stack with setup PSI0/PSIM residency disabled for copy-boundary diagnostics."
      ;;
    gpu_resident_stack_setup_psim_host)
      echo "Focused residency stack with only setup PSIM residency disabled for copy-boundary diagnostics."
      ;;
    gpu_resident_stack_psim_phase)
      echo "Focused residency stack plus cross-phase PSIM propagation residency."
      ;;
    gpu_resident_stack_hpsi_prop_psim_phase)
      echo "Focused residency stack plus HPSI-to-propagate and cross-phase PSIM propagation residency."
      ;;
    gpu_resident_stack_hpsi_prop_psim_switch)
      echo "Focused residency stack plus HPSI-to-propagate and PSIM-to-next-PSI0 switch residency."
      ;;
    gpu_resident_stack_cufft)
      echo "Focused residency stack plus threshold-gated native cuFFT for LIB\$FFTC8 calls."
      ;;
    gpu_resident_stack_cufft_force)
      echo "Focused residency stack plus force-all native cuFFT for LIB\$FFTC8 calls."
      ;;
    gpu_resident_stack_serial3dfft)
      echo "Focused residency stack plus opt-in single-rank full-grid 3D cuFFT path for PLANEWAVE\$FFT."
      ;;
    gpu_resident_stack_serial3dfft_force_dedpro)
      echo "Focused residency stack plus single-rank 3D cuFFT and opt-in force DEDPRO/PROFORCE device path."
      ;;
    gpu_resident_stack_serial3dfft_accmap)
      echo "Focused residency stack plus single-rank 3D cuFFT with device-side sparse/full-grid mapping."
      ;;
    gpu_resident_stack_serial3dfft_accmap_cache)
      echo "Focused residency stack plus cached device-side sparse/full-grid mapping for single-rank 3D cuFFT."
      ;;
    gpu_resident_stack_serial3dfft_accmap_vpsi_internal)
      echo "Device-side 3D FFT mapping diagnostic with resident VPSI real-space scratch."
      ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog)
      echo "Device-side 3D FFT mapping diagnostic that also keeps the HPSI RTOG output resident for the following HPSI consumers."
      ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal)
      echo "Device-side 3D FFT mapping diagnostic with HPSI RTOG output residency and resident VPSI real-space scratch."
      ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache)
      echo "Device-side 3D FFT mapping diagnostic with resident VPSI scratch and cached ACCMAP work arrays."
      ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache)
      echo "Device-side 3D FFT mapping diagnostic with cached ACCMAP work arrays and opt-in GPU density accumulation."
      ;;
    gpu_resident_stack_density_1cov_batch)
      echo "Current density-resident stack plus opt-in batched orthogonalization one-center overlap."
      ;;
    gpu_resident_stack_density_1cov_addoproj_cusolver_gram)
      echo "Current density-resident stack plus ADDOPROJ slice GEMMs and opt-in large-matrix cuSOLVER Gram-Cholesky."
      ;;
    gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj)
      echo "Current density-resident stack plus Gram-Cholesky and force-side resident ADDOPROJ cuBLAS."
      ;;
    gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum)
      echo "Full residency diagnostic combining HPSI, OPSI, DENMAT energy, persistent THIS%PROJ, and off-site device-pack accumulation."
      ;;
    gpu_resident_1coverlap_host)
      echo "Residency diagnostic that disables only the one-center overlap cuBLAS path."
      ;;
    *_invbatch_off)
      echo "cuBLAS diagnostic that disables batched inversion-symmetry scalarproducts."
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
    cusolver_standard)
      echo "cuSOLVER diagnostic: enables only standard DSYEVD/ZHEEVD offload."
      ;;
    cusolver_generalized*)
      echo "cuSOLVER diagnostic: enables only generalized DSYGVD/ZHEGVD offload."
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
    mpi_run=$(case_mpirun "${case_name}" "${exe}")
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
  if [[ -n "${CPPAW_CUSOLVER_ACC_STANDARD_MIN_N:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_ACC_STANDARD_MIN_N}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK=${CPPAW_CUSOLVER_ACC_CHECK}"
  fi
  if [[ -n "${CPPAW_CUSOLVER_ACC_CHECK_TOL:-}" ]]; then
    env_line="${env_line} CPPAW_CUSOLVER_ACC_CHECK_TOL=${CPPAW_CUSOLVER_ACC_CHECK_TOL}"
  fi
  echo "${env_line}"
}

cusolver_standard_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_FORCE_MIN_N:-1}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_OFF_MIN_N:-1000000000}"
  echo "${env_line}"
}

cusolver_generalized_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_OFF_MIN_N:-1000000000}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_FORCE_MIN_N:-1}"
  echo "${env_line}"
}

cusolver_generalized_conservative_env() {
  local env_line
  env_line=$(cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}")
  env_line="${env_line} CPPAW_CUSOLVER_ACC_STANDARD_MIN_N=${CPPAW_CUSOLVER_STANDARD_OFF_MIN_N:-1000000000}"
  env_line="${env_line} CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=${CPPAW_CUSOLVER_GENERALIZED_CONSERVATIVE_MIN_N:-256}"
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

inherited_accel_env() {
  local name value env_line=""
  while IFS='=' read -r name value; do
    case "${name}" in
      CPPAW_GPU_*|CPPAW_CUBLAS_ACC_*|CPPAW_CUSOLVER_ACC_*|CPPAW_CUFFT_ACC*|CPPAW_GRAM_CHOLESKY)
        if [[ "${value}" =~ ^[A-Za-z0-9_./:+-]+$ ]]; then
          env_line="${env_line:+${env_line} }${name}=${value}"
        fi
        ;;
    esac
  done < <(env)
  echo "${env_line}"
}

case_env() {
  case "$1" in
    cublas) cublas_env ;;
    cublas_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    cublas_invbatch_off) echo "$(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0" ;;
    cublas_conservative) cublas_conservative_env ;;
    cublas_projection_conservative) cublas_projection_conservative_env ;;
    cublas_overlap_conservative) cublas_overlap_conservative_env ;;
    cublas_addproduct_conservative) cublas_addproduct_conservative_env ;;
    cublas_matmul_conservative) cublas_matmul_conservative_env ;;
    cublas_off) echo "CPPAW_CUBLAS_ACC=0" ;;
    cusolver) cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}" ;;
    cusolver_standard) cusolver_standard_env ;;
    cusolver_generalized) cusolver_generalized_env ;;
    cusolver_generalized_conservative) cusolver_generalized_conservative_env ;;
    cusolver_conservative) cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}" ;;
    cusolver_off) echo "CPPAW_CUSOLVER_ACC=0" ;;
    cufft) cufft_env ;;
    cufft_force_all) cufft_force_env ;;
    cufft_off) echo "CPPAW_CUFFT_ACC=0" ;;
    gpu) cublas_env ;;
    gpu_nosync) echo "$(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_invbatch_off) echo "$(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0" ;;
    gpu_projection_conservative) cublas_projection_conservative_env ;;
    gpu_overlap_conservative) cublas_overlap_conservative_env ;;
    gpu_addproduct_conservative) cublas_addproduct_conservative_env ;;
    gpu_matmul_conservative) cublas_matmul_conservative_env ;;
    gpu_force_all) echo "$(cufft_force_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_conservative) echo "$(cublas_conservative_env) $(cusolver_env "${CPPAW_CUSOLVER_CONSERVATIVE_MIN_N:-256}")" ;;
    gpu_all) echo "$(cufft_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_nosync) echo "$(cufft_env) $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_invbatch_off) echo "$(cufft_env) $(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0 $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_3dfft) echo "$(cufft3d_env) $(cublas_env) $(cusolver_env "${CPPAW_CUSOLVER_ACC_MIN_N:-1}")" ;;
    gpu_all_off) echo "CPPAW_CUFFT_ACC=0 CPPAW_CUBLAS_ACC=0 CPPAW_CUSOLVER_ACC=0" ;;
    gpu_resident) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_nosync) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_resident_invbatch_off) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_INVERSION_BATCH=0" ;;
    gpu_resident_no_cusolver) echo "CPPAW_GPU_RESIDENCY=1 $(cublas_env) CPPAW_CUSOLVER_ACC=0" ;;
    gpu_resident_1coverlap) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_1COVERLAP=1 $(cublas_env)" ;;
    gpu_resident_1coverlap_batch) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_1COVERLAP=1 CPPAW_GPU_1COVERLAP_BATCH=1 $(cublas_env)" ;;
    gpu_resident_addoproj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_ADDOPROJ=1 CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO:-64} $(cublas_env)" ;;
    gpu_resident_pro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PRO_EXPANSION=0 $(cublas_env)" ;;
    gpu_resident_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_addpro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE=0 $(cublas_env)" ;;
    gpu_resident_addpro_hpsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_HPSI=0 $(cublas_env)" ;;
    gpu_resident_addpro_opsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_OPSI=0 $(cublas_env)" ;;
    gpu_resident_opsi_addpro_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_ADDPRO_CACHE_OPSI=0 $(cublas_env)" ;;
    gpu_resident_forcepsi_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_FORCE_PSI_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_orthoconst) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_CONST_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_orthox) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_orthox_off) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_orthox_nosync) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_ORTHO_X_RESIDENCY=1 $(cublas_env) CPPAW_CUBLAS_ACC_SYNC=0" ;;
    gpu_resident_opsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_denmat_energy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_denmat_energy_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_batch) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_batch) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_resident_hpsi_offden_blas) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 $(cublas_env)" ;;
    gpu_psim_propagate|gpu_resident_psim) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 $(cublas_env)" ;;
    gpu_hpsi_psim_propagate|gpu_resident_hpsi_psim) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 $(cublas_env)" ;;
    gpu_resident_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_psim_phase) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_proj) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_stack) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cublas_env)" ;;
    gpu_resident_stack_force_dedpro) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_psi0_ortho_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_ORTHO_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_psi0_prinfo_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_PRINFO_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_setup_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_SETUP_PSI_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_setup_psim_host) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_SETUP_PSIM_RESIDENCY=0 $(cublas_env)" ;;
    gpu_resident_stack_psim_phase) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_hpsi_prop_psim_phase) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_hpsi_prop_psim_switch) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1 CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1 $(cublas_env)" ;;
    gpu_resident_stack_cufft) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cufft_env) $(cublas_env)" ;;
    gpu_resident_stack_cufft_force) echo "CPPAW_GPU_RESIDENCY_STACK=1 $(cufft_force_env) $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_force_dedpro) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_cache) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_vpsi_internal) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_density_1cov_batch) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_GPU_1COVERLAP_BATCH=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_density_1cov_addoproj) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_GPU_1COVERLAP_BATCH=1 CPPAW_GPU_ORTHO_ADDOPROJ=1 CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO:-64} CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_density_1cov_addoproj_cusolver_gram) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_GPU_1COVERLAP_BATCH=1 CPPAW_GPU_ORTHO_ADDOPROJ=1 CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO:-64} CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY=1 CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N=${CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N:-4096} CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj) echo "CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_GPU_1COVERLAP_BATCH=1 CPPAW_GPU_ORTHO_ADDOPROJ=1 CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO:-64} CPPAW_GPU_FORCE_ADDOPROJ=1 CPPAW_GPU_FORCE_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_FORCE_ADDOPROJ_MIN_NPRO:-1} CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY=1 CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N=${CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N:-4096} CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack) echo "CPPAW_GPU_PROJECTION_STACK=1 CPPAW_GPU_ADDPRO_STACK=1 CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_FFT_SERIAL_3D_ACC_CACHE=1 CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1 CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1 CPPAW_GPU_1COVERLAP_BATCH=1 CPPAW_GPU_ORTHO_ADDOPROJ=1 CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO:-64} CPPAW_GPU_FORCE_ADDOPROJ=1 CPPAW_GPU_FORCE_ADDOPROJ_MIN_NPRO=${CPPAW_GPU_FORCE_ADDOPROJ_MIN_NPRO:-1} CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY=1 CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N=${CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N:-4096} CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=${CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS:-0} $(cublas_env)" ;;
    gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_HPSI_RESIDENCY=1 CPPAW_GPU_OPSI_RESIDENCY=1 CPPAW_GPU_PROJ_RESIDENCY=1 CPPAW_GPU_DENMAT_ENERGY=1 CPPAW_GPU_DENMAT_MINFLOP=${CPPAW_GPU_DENMAT_MINFLOP:-1} CPPAW_GPU_OFFDEN_LOCAL=1 CPPAW_GPU_OFFDEN_CUBLAS=1 CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1 CPPAW_GPU_OFFDEN_DEVICE_PACK=1 CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1 CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=${CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP:-1} CPPAW_GPU_OFFDEN_BATCH_SIZE=${CPPAW_GPU_OFFDEN_BATCH_SIZE:-64} $(cublas_env)" ;;
    gpu_resident_1coverlap_host) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GPU_1COVERLAP=0 $(cublas_env)" ;;
    gpu_resident_gram_cholesky) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GRAM_CHOLESKY=1 $(cublas_env)" ;;
    gpu_resident_gram_legacy) echo "CPPAW_GPU_RESIDENCY=1 CPPAW_GRAM_CHOLESKY=0 $(cublas_env)" ;;
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

case_runtime_env() {
  local case_name=$1
  local exe=${2:-}
  local libdir
  case "${case_name}" in
    cpu)
      libdir=$(cpu_mpi_libdir_for_exe "${exe}")
      if [[ "${RANKS}" -gt 1 && -n "${libdir}" ]]; then
        echo "LD_LIBRARY_PATH=${libdir}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
      fi
      ;;
    *) echo "" ;;
  esac
}

combined_env() {
  local env_line=""
  local part
  for part in "$@"; do
    [[ -n "${part}" ]] || continue
    env_line="${env_line:+${env_line} }${part}"
  done
  echo "${env_line}"
}

iso_now() {
  date -u +%Y-%m-%dT%H:%M:%SZ
}

prepare_case() {
  local dir=$1
  local skala_device skala_check
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
  if grep -q '@SKALA_MODEL@' "${dir}/${TEST}.cntl"; then
    if [[ -z "${SKALA_MODEL:-}" ]]; then
      echo "SKALA_MODEL is required by ${CNTL_FILE}" >&2
      return 1
    fi
    skala_device=${SKALA_DEVICE:-AUTO}
    skala_check=${SKALA_CHECK:-F}
    case "${skala_check}" in
      T|F) ;;
      *)
        echo "SKALA_CHECK must be T or F" >&2
        return 1
        ;;
    esac
    SKALA_MODEL="${SKALA_MODEL}" SKALA_DEVICE="${skala_device}" \
      SKALA_CHECK="${skala_check}" perl -0pi -e \
      's/\@SKALA_MODEL\@/$ENV{SKALA_MODEL}/g;
       s/\@SKALA_DEVICE\@/$ENV{SKALA_DEVICE}/g;
       s/\@SKALA_CHECK\@/$ENV{SKALA_CHECK}/g' "${dir}/${TEST}.cntl"
  fi
  if [[ -n "${EMPTY_BANDS:-}" ]]; then
    perl -0pi -e "s/EMPTY\\s*=\\s*\\d+/EMPTY=${EMPTY_BANDS}/" "${dir}/${TEST}.strc"
  fi
}

capture_metadata() {
  local missing_cases=0
  {
    echo "date=$(iso_now)"
    echo "hostname=$(hostname)"
    echo "test=${TEST}"
    echo "empty_bands=${EMPTY_BANDS:-}"
    echo "nsteps=${NSTEPS}"
    echo "ranks=${RANKS}"
    echo "cases=${CASES}"
    if [[ -n "${EXPECTED_ENERGY}" ]]; then
      echo "expected_energy=${EXPECTED_ENERGY}"
      echo "energy_tol=${ENERGY_TOL}"
    fi
    echo "threads=OMP_NUM_THREADS=${OMP_NUM_THREADS} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS} BLIS_NUM_THREADS=${BLIS_NUM_THREADS} VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS} NVPL_NUM_THREADS=${NVPL_NUM_THREADS}"
    [[ -n "${PKG_CONFIG_PATH:-}" ]] && echo "pkg_config_path=${PKG_CONFIG_PATH}"
    [[ -n "${LD_LIBRARY_PATH:-}" ]] && echo "ld_library_path=${LD_LIBRARY_PATH}"
    echo
    uname -a
    echo
    for case_name in ${CASES}; do
      exe=$(if [[ "${RANKS}" -gt 1 ]]; then parallel_exe "${case_name}"; else serial_exe "${case_name}"; fi)
      echo "case=${case_name}"
      echo "exe=${exe}"
      if [[ "${RANKS}" -gt 1 ]]; then
        echo "mpirun=$(case_mpirun "${case_name}" "${exe}")"
      fi
      note=$(case_note "${case_name}")
      [[ -n "${note}" ]] && echo "note=${note}"
      extra_env=$(inherited_accel_env)
      [[ -n "${extra_env}" ]] && echo "inherited_accel_env=${extra_env}"
      runtime_env=$(case_runtime_env "${case_name}" "${exe}")
      [[ -n "${runtime_env}" ]] && echo "runtime_env=${runtime_env}"
      case_env_line=$(case_env "${case_name}")
      [[ -n "${case_env_line}" ]] && echo "case_env=${case_env_line}"
      env_line=$(combined_env "${extra_env}" "${runtime_env}" "${case_env_line}")
      [[ -n "${env_line}" ]] && echo "planned_env=${env_line}"
      cmd=$(run_command "${case_name}" "${exe}")
      echo "planned_command=${cmd}"
      time_part="${TIME_CMD:+${TIME_CMD} -p }"
      timeout_part="${TIMEOUT_PREFIX:+${TIMEOUT_PREFIX} }"
      if [[ -n "${env_line}" ]]; then
        echo "planned_full_command=${time_part}env CPPAW_ACCEL_PROFILE_FILE=${case_name}_profile ${env_line} ${timeout_part}${cmd}"
      else
        echo "planned_full_command=${time_part}env CPPAW_ACCEL_PROFILE_FILE=${case_name}_profile ${timeout_part}${cmd}"
      fi
      if [[ -x "${exe}" ]]; then
        if [[ -n "${runtime_env}" ]]; then
          # shellcheck disable=SC2086
          env ${runtime_env} ldd "${exe}" 2>/dev/null | grep -E "blas|lapack|fftw|cufft|cusolver|cublas|nvpl|openblas|libmpi" || true
        else
          ldd "${exe}" 2>/dev/null | grep -E "blas|lapack|fftw|cufft|cusolver|cublas|nvpl|openblas|libmpi" || true
        fi
      else
        echo "missing"
        missing_cases=$((missing_cases+1))
      fi
      echo
    done
  } > "${RUN_ROOT}/metadata.txt" 2>&1
  if [[ "${missing_cases}" -gt 0 \
      && ( "${REQUIRE_CASES}" == "yes" \
           || "${REQUIRE_CASES}" == "true" \
           || "${REQUIRE_CASES}" == "1" ) ]]; then
    echo "Required benchmark cases missing executable: ${missing_cases}" >&2
    echo "See ${RUN_ROOT}/metadata.txt for case details." >&2
    return 1
  fi
}

mkdir -p "${RUN_ROOT}" "${HERE}/runs"
echo "${RUN_ROOT}" > "${HERE}/runs/latest"
capture_metadata
case "${DRY_RUN}" in
  yes|true|1)
    echo "Benchmark dry-run metadata: ${RUN_ROOT}"
    exit 0
    ;;
esac

for case_name in ${CASES}; do
  exe=$(if [[ "${RANKS}" -gt 1 ]]; then parallel_exe "${case_name}"; else serial_exe "${case_name}"; fi)
  if [[ ! -x "${exe}" ]]; then
    if [[ "${REQUIRE_CASES}" == "yes" || "${REQUIRE_CASES}" == "true" || "${REQUIRE_CASES}" == "1" ]]; then
      echo "Required case missing executable: ${case_name}" >&2
      echo "  exe=${exe}" >&2
      exit 1
    fi
    echo "Skipping ${case_name}: executable not found: ${exe}" >&2
    continue
  fi
  for repeat in $(seq 1 "${REPEATS}"); do
    run_dir="${RUN_ROOT}/${case_name}/rep$(printf "%02d" "${repeat}")"
    prepare_case "${run_dir}"
    (
      cd "${run_dir}"
      env_line=$(combined_env "$(inherited_accel_env)" \
        "$(case_runtime_env "${case_name}" "${exe}")" "$(case_env "${case_name}")")
      {
        echo "case=${case_name}"
        echo "repeat=${repeat}"
        echo "nsteps=${NSTEPS}"
        echo "ranks=${RANKS}"
        echo "exe=${exe}"
        echo "env=${env_line}"
        if [[ -n "${EXPECTED_ENERGY}" ]]; then
          echo "expected_energy=${EXPECTED_ENERGY}"
          echo "energy_tol=${ENERGY_TOL}"
        fi
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
