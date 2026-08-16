#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)

cd "${ROOT}"

bash -n paw_install
bash -n src/Buildtools/defaultparmfile
bash -n src/Buildtools/paw_build.sh
bash -n src/Buildtools/paw_fcflags.sh
bash -n src/Buildtools/paw_srclist.sh
bash -n tests/profile/si64/run_benchmark.sh
bash -n tests/profile/si64/case_recommendations.sh
bash -n tests/profile/si64/run_cusolver_focus.sh
bash -n tests/profile/si64/run_followup.sh
bash -n tests/profile/si64/run_gap_profile_night.sh
bash -n tests/profile/si64/run_gpu_resource_comparison.sh
bash -n tests/profile/si64/run_large_bands.sh
bash -n tests/profile/si64/run_large_bands_long.sh
bash -n tests/profile/si64/run_gpu_exploration.sh
bash -n tests/profile/si64/run_gpu_library_matrix.sh
bash -n tests/profile/si64/run_nvhpc_standard.sh
bash -n tests/profile/si64/run_nsys.sh
bash -n tests/profile/si64/run_offden_focus.sh
bash -n tests/profile/si64/run_overnight.sh
bash -n tests/profile/si64/run_psim_focus.sh
bash -n tests/profile/si64/run_psim_lifecycle.sh
bash -n tests/profile/si64/run_vpsi_boundary.sh
bash -n src/Tools/Scripts/paw_cuda_aware_mpi_probe.sh
bash -n src/Tools/Scripts/paw_gpu_capabilities.sh

python3 -m py_compile \
  tests/profile/si64/profile_summary.py \
  tests/profile/si64/profile_copy_rows.py \
  tests/profile/si64/check_benchmark_tools.py \
  tests/profile/si64/benchmark_compare.py \
  tests/profile/si64/benchmark_summary.py \
  tests/profile/si64/benchmark_markdown.py \
  tests/profile/si64/nsys_sql_summary.py

python3 tests/profile/si64/check_benchmark_tools.py

tmpdir=$(mktemp -d)
trap 'rm -rf "${tmpdir}"' EXIT

src/Tools/Scripts/paw_gpu_capabilities.sh > "${tmpdir}/gpu_capabilities.txt"
grep -q "^host_fftw=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^host_fftw_ld_library_path=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^host_fftw_pkg_config_path=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^host_blas_lapack=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_cpu_cases=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_gpu_cases=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_gpu_diagnostic_cases=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_large_band_gpu_cases=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_resource_cases=" "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_standard_command=cd tests/profile/si64" \
  "${tmpdir}/gpu_capabilities.txt"
grep -q "^recommended_large_band_command=cd tests/profile/si64" \
  "${tmpdir}/gpu_capabilities.txt"

cat > "${tmpdir}/fake_gpu_capabilities.txt" <<'EOF'
host_fftw=yes path=/tmp/fftw/lib/libfftw3 include=/tmp/fftw/include/fftw3.f03
host_fftw_ld_library_path=/tmp/fftw/lib
host_fftw_pkg_config_path=/tmp/fftw/lib/pkgconfig
recommended_cpu_cases=cpu nvhpc_cpu
recommended_gpu_cases=gpu_resident_stack gpu_resident_off gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_cufft
recommended_gpu_diagnostic_cases=gpu_resident_stack_force_dedpro gpu_resident_nosync gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache
recommended_large_band_gpu_cases=gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_force_dedpro gpu_resident_stack_serial3dfft_accmap_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache
recommended_resource_cases=cpu nvhpc_cpu gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack
EOF

cat > "${tmpdir}/fake_nohost_capabilities.txt" <<'EOF'
host_fftw=no
host_fftw_ld_library_path=none
host_fftw_pkg_config_path=none
host_blas_lapack=yes path=/opt/nvidia/hpc_sdk/Linux_x86_64/2024/compilers/lib/lib{blas,lapack}.so
recommended_cpu_reason=no_host_fftw_runtime_found
recommended_gpu_reason=no_host_fftw_runtime_found
recommended_cpu_cases=none
recommended_gpu_cases=none
recommended_gpu_diagnostic_cases=none
recommended_large_band_gpu_cases=none
recommended_resource_cases=none
EOF

PKG_CONFIG_PATH= \
LD_LIBRARY_PATH= \
CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_gpu_capabilities.txt" \
bash -c '. tests/profile/si64/case_recommendations.sh
cppaw_apply_capability_env
case ":${PKG_CONFIG_PATH}:" in *":/tmp/fftw/lib/pkgconfig:"*) ;; *) exit 1 ;; esac
case ":${LD_LIBRARY_PATH}:" in *":/tmp/fftw/lib:"*) ;; *) exit 1 ;; esac'

DRY_RUN=yes \
  RUN_ROOT="${tmpdir}/dry-run" \
  CASES="cufft_force_all gpu_all_resident_stack gpu_all_resident_stack_cufft gpu_resident_stack gpu_resident_stack_force_dedpro gpu_resident_stack_psi0_ortho_host gpu_resident_stack_psi0_prinfo_host gpu_resident_stack_setup_psim_host gpu_resident_stack_hpsi_prop_psim_phase gpu_resident_stack_hpsi_prop_psim_switch gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_force_dedpro gpu_resident_stack_serial3dfft_accmap gpu_resident_stack_serial3dfft_accmap_cache gpu_resident_stack_serial3dfft_accmap_vpsi_internal gpu_resident_stack_serial3dfft_accmap_hpsi_rtog gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_no_cufft" \
  tests/profile/si64/run_benchmark.sh > "${tmpdir}/dry-run.out"
grep -q "Benchmark dry-run metadata:" "${tmpdir}/dry-run.out"
grep -q "case=cufft_force_all" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_all_resident_stack" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_all_resident_stack_cufft" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack" "${tmpdir}/dry-run/metadata.txt"
grep -q "^empty_bands=$" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_force_dedpro" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_psi0_ortho_host" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_ORTHO_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_PSI0_ORTHO_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_psi0_prinfo_host" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_PRINFO_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_PSI0_PRINFO_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_setup_psim_host" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_SETUP_PSIM_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_SETUP_PSIM_RESIDENCY=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_hpsi_prop_psim_phase" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_hpsi_prop_psim_switch" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSIM_PROPAGATE=1 CPPAW_GPU_PSIM_PHASE_RESIDENCY=1 CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1 CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_force_dedpro" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1.*CPPAW_FFT_SERIAL_3D=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D_ACC_MAP=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_cache" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D_ACC_CACHE=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_vpsi_internal" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_FFT_SERIAL_3D_ACC_MAP=1 CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D_ACC_CACHE=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_FORCE_ADDOPROJ=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_FORCE_ADDOPROJ_MIN_NPRO=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_PROJECTION_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_ADDPRO_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_FFT_SERIAL_3D=1 CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_no_cufft" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_CUFFT_ACC=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_CUFFT_ACC=0" "${tmpdir}/dry-run/metadata.txt"

DRY_RUN=yes \
  RUN_ROOT="${tmpdir}/dry-run-nsteps3" \
  NSTEPS=3 \
  CASES="gpu_resident_stack" \
  tests/profile/si64/run_benchmark.sh > "${tmpdir}/dry-run-nsteps3.out"
! grep -q "^expected_energy=" "${tmpdir}/dry-run-nsteps3/metadata.txt"

DRY_RUN=yes \
  RUN_ROOT="${tmpdir}/dry-run-nsteps3-explicit" \
  NSTEPS=3 \
  EXPECTED_ENERGY=208.886424 \
  CASES="gpu_resident_stack" \
  tests/profile/si64/run_benchmark.sh \
    > "${tmpdir}/dry-run-nsteps3-explicit.out"
grep -q "expected_energy=208.886424" \
  "${tmpdir}/dry-run-nsteps3-explicit/metadata.txt"

set +e
DRY_RUN=yes \
  REQUIRE_CASES=yes \
  RUN_ROOT="${tmpdir}/dry-run-required" \
  CASES="gpu_resident_stack" \
  tests/profile/si64/run_benchmark.sh > "${tmpdir}/dry-run-required.out" 2>&1
dry_run_required_status=$?
set -e
if grep -q "^missing$" "${tmpdir}/dry-run-required/metadata.txt"; then
  test "${dry_run_required_status}" -ne 0
  grep -q "Required benchmark cases missing executable:" \
    "${tmpdir}/dry-run-required.out"
else
  test "${dry_run_required_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  NVHPC_STANDARD_ROOT="${tmpdir}/standard-dry-run" \
  GPU_CASES="gpu_resident_stack" \
  CPU_CASES="" \
  tests/profile/si64/run_nvhpc_standard.sh > "${tmpdir}/standard-dry-run.out" 2>&1
standard_dry_run_status=$?
set -e
grep -q "DRY-RUN skip pre-build binary check" \
  "${tmpdir}/standard-dry-run/nvhpc_standard.log"
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" \
  "${tmpdir}/standard-dry-run/gpu_1rank/metadata.txt"
if grep -q "^missing$" "${tmpdir}/standard-dry-run/gpu_1rank/metadata.txt"; then
  test "${standard_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/standard-dry-run/nvhpc_standard.log"
else
  test "${standard_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_gpu_capabilities.txt" \
  NVHPC_STANDARD_ROOT="${tmpdir}/standard-auto-dry-run" \
  GPU_CASES=auto \
  CPU_CASES=auto \
  tests/profile/si64/run_nvhpc_standard.sh \
    > "${tmpdir}/standard-auto-dry-run.out" 2>&1
standard_auto_dry_run_status=$?
set -e
grep -q "selected_cases gpu='gpu_resident_stack gpu_resident_off gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_cufft' cpu='cpu nvhpc_cpu'" \
  "${tmpdir}/standard-auto-dry-run/nvhpc_standard.log"
grep -q "^recommended_gpu_cases=gpu_resident_stack gpu_resident_off gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_cufft" \
  "${tmpdir}/standard-auto-dry-run/gpu_capabilities.txt"
grep -q "case=gpu_resident_stack" \
  "${tmpdir}/standard-auto-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_off" \
  "${tmpdir}/standard-auto-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack" \
  "${tmpdir}/standard-auto-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_cufft" \
  "${tmpdir}/standard-auto-dry-run/gpu_1rank/metadata.txt"
grep -q "case=cpu" "${tmpdir}/standard-auto-dry-run/cpu_1rank/metadata.txt"
grep -q "case=nvhpc_cpu" \
  "${tmpdir}/standard-auto-dry-run/cpu_8rank_ref/metadata.txt"
if grep -q "^missing$" "${tmpdir}/standard-auto-dry-run/gpu_1rank/metadata.txt" \
    || grep -q "^missing$" "${tmpdir}/standard-auto-dry-run/cpu_1rank/metadata.txt" \
    || grep -q "^missing$" "${tmpdir}/standard-auto-dry-run/cpu_8rank_ref/metadata.txt"; then
  test "${standard_auto_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/standard-auto-dry-run/nvhpc_standard.log"
else
  test "${standard_auto_dry_run_status}" -eq 0
fi

DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_nohost_capabilities.txt" \
  NVHPC_STANDARD_ROOT="${tmpdir}/standard-nohost-dry-run" \
  GPU_CASES=auto \
  CPU_CASES=auto \
  tests/profile/si64/run_nvhpc_standard.sh \
    > "${tmpdir}/standard-nohost-dry-run.out" 2>&1
grep -q "selected_cases gpu='none' cpu='none'" \
  "${tmpdir}/standard-nohost-dry-run/nvhpc_standard.log"
grep -q "SKIP all suites empty case lists" \
  "${tmpdir}/standard-nohost-dry-run/nvhpc_standard.log"
test ! -e "${tmpdir}/standard-nohost-dry-run/gpu_1rank/metadata.txt"
test ! -e "${tmpdir}/standard-nohost-dry-run/cpu_1rank/metadata.txt"

DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_gpu_capabilities.txt" \
  OVERNIGHT_ROOT="${tmpdir}/overnight-auto-dry-run" \
  LONG_NSTEPS=1 \
  SCALING_NSTEPS=1 \
  THRESHOLD_NSTEPS=1 \
  LONG_REPEATS=1 \
  SCALING_REPEATS=1 \
  THRESHOLD_REPEATS=1 \
  tests/profile/si64/run_overnight.sh \
    > "${tmpdir}/overnight-auto-dry-run.out" 2>&1
grep -q "selected_cases main='cpu nvhpc_cpu gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack' scaling='cpu nvhpc_cpu gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack'.*threshold='gpu_resident_stack' nsys='gpu_resident_stack'" \
  "${tmpdir}/overnight-auto-dry-run/overnight.log"
grep -q "case=gpu_resident_stack" \
  "${tmpdir}/overnight-auto-dry-run/main_1steps_4ranks/metadata.txt"
grep -q "case=cpu" \
  "${tmpdir}/overnight-auto-dry-run/scaling_1steps_1ranks/metadata.txt"
grep -q "case=nvhpc_cpu" \
  "${tmpdir}/overnight-auto-dry-run/scaling_1steps_4ranks/metadata.txt"
grep -q "case=gpu_resident_stack" \
  "${tmpdir}/overnight-auto-dry-run/threshold_1e7_1steps_4ranks/metadata.txt"
grep -q "DRY-RUN skip suite=nsys_nstep5_1ranks case=gpu_resident_stack" \
  "${tmpdir}/overnight-auto-dry-run/overnight.log"
test "$(cat "${tmpdir}/overnight-auto-dry-run/nsys_nstep5_1ranks.status")" = "dry-run"

DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_nohost_capabilities.txt" \
  OVERNIGHT_ROOT="${tmpdir}/overnight-nohost-dry-run" \
  LONG_NSTEPS=1 \
  SCALING_NSTEPS=1 \
  THRESHOLD_NSTEPS=1 \
  LONG_REPEATS=1 \
  SCALING_REPEATS=1 \
  THRESHOLD_REPEATS=1 \
  tests/profile/si64/run_overnight.sh \
    > "${tmpdir}/overnight-nohost-dry-run.out" 2>&1
grep -q "selected_cases main='none' scaling='none'.*threshold='none' nsys='none'" \
  "${tmpdir}/overnight-nohost-dry-run/overnight.log"
grep -q "SKIP  suite=main_1steps_4ranks empty case list" \
  "${tmpdir}/overnight-nohost-dry-run/overnight.log"
grep -q "SKIP  suite=threshold empty case list" \
  "${tmpdir}/overnight-nohost-dry-run/overnight.log"
grep -q "SKIP  suite=nsys disabled_or_empty_case" \
  "${tmpdir}/overnight-nohost-dry-run/overnight.log"
test ! -e "${tmpdir}/overnight-nohost-dry-run/main_1steps_4ranks/metadata.txt"

set +e
DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_gpu_capabilities.txt" \
  GPU_RESOURCE_COMPARISON_ROOT="${tmpdir}/resource-dry-run" \
  GPU_CASES=auto \
  CPU_CASES="" \
  tests/profile/si64/run_gpu_resource_comparison.sh \
    > "${tmpdir}/resource-dry-run.out" 2>&1
resource_dry_run_status=$?
set -e
grep -q "selected_cases gpu='gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_force_dedpro gpu_resident_stack_serial3dfft_accmap_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache' cpu='none'" \
  "${tmpdir}/resource-dry-run/nvhpc_standard.log"
grep -q "^recommended_large_band_gpu_cases=gpu_resident_stack gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_serial3dfft gpu_resident_stack_serial3dfft_force_dedpro gpu_resident_stack_serial3dfft_accmap_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache" \
  "${tmpdir}/resource-dry-run/gpu_capabilities.txt"
grep -q "case=gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_force_dedpro" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "^empty_bands=2048$" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "^pkg_config_path=/tmp/fftw/lib/pkgconfig" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "^ld_library_path=/tmp/fftw/lib" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1.*CPPAW_FFT_SERIAL_3D=1" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_cache" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D_ACC_CACHE=1" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "planned_full_command=.*CPPAW_FFT_SERIAL_3D_ACC_CACHE=1" \
  "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"
grep -q "SKIP  suite=cpu_1rank empty case list" \
  "${tmpdir}/resource-dry-run/nvhpc_standard.log"
test ! -e "${tmpdir}/resource-dry-run/cpu_1rank/metadata.txt"
if grep -q "^missing$" "${tmpdir}/resource-dry-run/gpu_1rank/metadata.txt"; then
  test "${resource_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/resource-dry-run/nvhpc_standard.log"
else
  test "${resource_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  GPU_LIBRARY_MATRIX_ROOT="${tmpdir}/matrix-dry-run" \
  GPU_CASES="gpu_resident_stack" \
  CPU_CASES="" \
  tests/profile/si64/run_gpu_library_matrix.sh > "${tmpdir}/matrix-dry-run.out" 2>&1
matrix_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" \
  "${tmpdir}/matrix-dry-run/gpu_1rank/metadata.txt"
grep -q "SKIP  suite=cpu_1rank empty case list" \
  "${tmpdir}/matrix-dry-run/nvhpc_standard.log"
test ! -e "${tmpdir}/matrix-dry-run/cpu_1rank/metadata.txt"
test ! -e "${tmpdir}/matrix-dry-run/cpu_8rank_ref/metadata.txt"
if grep -q "^missing$" "${tmpdir}/matrix-dry-run/gpu_1rank/metadata.txt"; then
  test "${matrix_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/matrix-dry-run/nvhpc_standard.log"
else
  test "${matrix_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  GPU_EXPLORATION_ROOT="${tmpdir}/exploration-dry-run" \
  GPU_CASES="gpu_resident_stack" \
  CPU_CASES="" \
  tests/profile/si64/run_gpu_exploration.sh > "${tmpdir}/exploration-dry-run.out" 2>&1
exploration_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" \
  "${tmpdir}/exploration-dry-run/one_rank_gpu/metadata.txt"
grep -q "SKIP  suite=one_rank_cpu empty case list" \
  "${tmpdir}/exploration-dry-run/gpu_exploration.log"
test ! -e "${tmpdir}/exploration-dry-run/one_rank_cpu/metadata.txt"
test ! -e "${tmpdir}/exploration-dry-run/eight_rank_cpu/metadata.txt"
if grep -q "^missing$" "${tmpdir}/exploration-dry-run/one_rank_gpu/metadata.txt"; then
  test "${exploration_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/exploration-dry-run/gpu_exploration.log"
else
  test "${exploration_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  CPPAW_GPU_CAPABILITIES_FILE="${tmpdir}/fake_gpu_capabilities.txt" \
  GPU_EXPLORATION_ROOT="${tmpdir}/exploration-auto-dry-run" \
  GPU_CASES=auto \
  CPU_CASES="" \
  tests/profile/si64/run_gpu_exploration.sh \
    > "${tmpdir}/exploration-auto-dry-run.out" 2>&1
exploration_auto_dry_run_status=$?
set -e
grep -q "selected_cases gpu='gpu_resident_stack gpu_resident_off gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack gpu_resident_stack_cufft gpu_resident_stack_force_dedpro gpu_resident_nosync gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache' cpu='none'" \
  "${tmpdir}/exploration-auto-dry-run/gpu_exploration.log"
grep -q "case=gpu_resident_stack_density_1cov_addoproj_cusolver_gram_force_addoproj_projaddpro_stack" \
  "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"
grep -q "case=gpu_resident_stack_cufft" \
  "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"
grep -q "case=gpu_resident_stack_force_dedpro" \
  "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"
grep -q "case=gpu_resident_nosync" \
  "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"
grep -q "case=gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache" \
  "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"
grep -q "SKIP  suite=one_rank_cpu empty case list" \
  "${tmpdir}/exploration-auto-dry-run/gpu_exploration.log"
if grep -q "^missing$" \
    "${tmpdir}/exploration-auto-dry-run/one_rank_gpu/metadata.txt"; then
  test "${exploration_auto_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/exploration-auto-dry-run/gpu_exploration.log"
else
  test "${exploration_auto_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  FOLLOWUP_ROOT="${tmpdir}/followup-dry-run" \
  NSTEPS_LIST=1 \
  EMPTY_BANDS_LIST=128 \
  GPU_CASES="gpu_resident_stack" \
  ONE_RANK_CPU_CASES="" \
  CPU_CASES="" \
  tests/profile/si64/run_followup.sh > "${tmpdir}/followup-dry-run.out" 2>&1
followup_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" \
  "${tmpdir}/followup-dry-run/empty128_nstep1_1rank_gpu/metadata.txt"
grep -q "SKIP  suite=empty128_nstep1_1rank_cpu empty case list" \
  "${tmpdir}/followup-dry-run/followup.log"
test ! -e "${tmpdir}/followup-dry-run/empty128_nstep1_1rank_cpu/metadata.txt"
test ! -e "${tmpdir}/followup-dry-run/empty128_nstep1_8rank_cpu_ref/metadata.txt"
if grep -q "^missing$" \
    "${tmpdir}/followup-dry-run/empty128_nstep1_1rank_gpu/metadata.txt"; then
  test "${followup_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/followup-dry-run/followup.log"
else
  test "${followup_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  CUSOLVER_FOCUS_ROOT="${tmpdir}/cusolver-dry-run" \
  EMPTY_BANDS_LIST=128 \
  CUSOLVER_CASES="cusolver_generalized" \
  CPU_CASES="" \
  tests/profile/si64/run_cusolver_focus.sh > "${tmpdir}/cusolver-dry-run.out" 2>&1
cusolver_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N=1" \
  "${tmpdir}/cusolver-dry-run/empty128_1rank_cusolver/metadata.txt"
grep -q "SKIP  suite=empty128_1rank_cpu empty case list" \
  "${tmpdir}/cusolver-dry-run/cusolver_focus.log"
test ! -e "${tmpdir}/cusolver-dry-run/empty128_1rank_cpu/metadata.txt"
test ! -e "${tmpdir}/cusolver-dry-run/empty128_8rank_cpu_ref/metadata.txt"
if grep -q "^missing$" \
    "${tmpdir}/cusolver-dry-run/empty128_1rank_cusolver/metadata.txt"; then
  test "${cusolver_dry_run_status}" -ne 0
  grep -q "FAILED suites=" "${tmpdir}/cusolver-dry-run/cusolver_focus.log"
else
  test "${cusolver_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  RUN_ROOT_BASE="${tmpdir}/offden-dry-run" \
  OFFDEN_CASES="gpu_resident_hpsi_offden_cublas_devicepack" \
  EMPTY_BANDS=128 \
  RUN_SHARED_GPU=no \
  tests/profile/si64/run_offden_focus.sh > "${tmpdir}/offden-dry-run.out" 2>&1
offden_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_OFFDEN_DEVICE_PACK=1" \
  "${tmpdir}/offden-dry-run-128-1r/metadata.txt"
if grep -q "^missing$" "${tmpdir}/offden-dry-run-128-1r/metadata.txt"; then
  test "${offden_dry_run_status}" -ne 0
else
  test "${offden_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  RUN_ROOT_BASE="${tmpdir}/psim-focus-dry-run" \
  PSIM_CASES="gpu_resident_hpsi_opsi_psim_phase" \
  EMPTY_BANDS=128 \
  RUN_SHARED_GPU=no \
  RUN_LARGE_GPU=no \
  tests/profile/si64/run_psim_focus.sh > "${tmpdir}/psim-focus-dry-run.out" 2>&1
psim_focus_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_PSIM_PHASE_RESIDENCY=1" \
  "${tmpdir}/psim-focus-dry-run-128-1r/metadata.txt"
if grep -q "^missing$" "${tmpdir}/psim-focus-dry-run-128-1r/metadata.txt"; then
  test "${psim_focus_dry_run_status}" -ne 0
else
  test "${psim_focus_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  RUN_ROOT_BASE="${tmpdir}/psim-lifecycle-dry-run" \
  PSIM_LIFECYCLE_CASES="gpu_resident_hpsi_opsi_psim_phase" \
  NSTEPS_LIST=1 \
  EMPTY_BANDS_LIST=128 \
  RUN_SHARED_GPU=no \
  RUN_LARGE_GPU=no \
  tests/profile/si64/run_psim_lifecycle.sh > "${tmpdir}/psim-lifecycle-dry-run.out" 2>&1
psim_lifecycle_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_GPU_PSIM_PHASE_RESIDENCY=1" \
  "${tmpdir}/psim-lifecycle-dry-run-empty128-nstep1-1r/metadata.txt"
if grep -q "^missing$" \
    "${tmpdir}/psim-lifecycle-dry-run-empty128-nstep1-1r/metadata.txt"; then
  test "${psim_lifecycle_dry_run_status}" -ne 0
else
  test "${psim_lifecycle_dry_run_status}" -eq 0
fi

set +e
DRY_RUN=yes \
  RUN_ROOT_BASE="${tmpdir}/vpsi-dry-run" \
  VPSI_BOUNDARY_CASES="gpu_resident_stack_cufft_force" \
  NSTEPS_LIST=1 \
  EMPTY_BANDS_LIST=128 \
  RUN_SHARED_GPU=no \
  RUN_LARGE_GPU=no \
  tests/profile/si64/run_vpsi_boundary.sh > "${tmpdir}/vpsi-dry-run.out" 2>&1
vpsi_dry_run_status=$?
set -e
grep -q "planned_full_command=.*CPPAW_CUFFT_ACC_MIN_ELEMENTS=0" \
  "${tmpdir}/vpsi-dry-run-empty128-nstep1-1r/metadata.txt"
if grep -q "^missing$" "${tmpdir}/vpsi-dry-run-empty128-nstep1-1r/metadata.txt"; then
  test "${vpsi_dry_run_status}" -ne 0
else
  test "${vpsi_dry_run_status}" -eq 0
fi

test -f tests/profile/si64/si64.cntl
test -f tests/profile/si64/si64.strc
test -f tests/profile/si64/si64_bands.cntl

git diff --check
