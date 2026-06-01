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
bash -n tests/profile/si64/run_cusolver_focus.sh
bash -n tests/profile/si64/run_followup.sh
bash -n tests/profile/si64/run_gap_profile_night.sh
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
DRY_RUN=yes \
  RUN_ROOT="${tmpdir}/dry-run" \
  CASES="gpu_resident_stack gpu_no_cufft" \
  tests/profile/si64/run_benchmark.sh > "${tmpdir}/dry-run.out"
grep -q "Benchmark dry-run metadata:" "${tmpdir}/dry-run.out"
grep -q "case=gpu_resident_stack" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_GPU_RESIDENCY_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_GPU_RESIDENCY_STACK=1" "${tmpdir}/dry-run/metadata.txt"
grep -q "case=gpu_no_cufft" "${tmpdir}/dry-run/metadata.txt"
grep -q "case_env=CPPAW_CUFFT_ACC=0" "${tmpdir}/dry-run/metadata.txt"
grep -q "planned_full_command=.*CPPAW_CUFFT_ACC=0" "${tmpdir}/dry-run/metadata.txt"

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

test -f tests/profile/si64/si64.cntl
test -f tests/profile/si64/si64.strc
test -f tests/profile/si64/si64_bands.cntl

git diff --check
