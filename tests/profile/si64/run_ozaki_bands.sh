#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

export TEST=${TEST:-si64_bands}
export FOLLOWUP_ROOT=${OZAKI_ROOT:-"${HERE}/runs/${TEST}-ozaki-$(date +%Y%m%d-%H%M%S)"}
export NSTEPS_LIST=${NSTEPS_LIST:-1}
export EMPTY_BANDS_LIST=${EMPTY_BANDS_LIST:-"1024 2048 4096"}
export REPEATS=${REPEATS:-1}
export GPU_RANKS=${GPU_RANKS:-1}
export CPU_RANKS=${CPU_RANKS:-8}
export CPPAW_CUBLAS_FP64_WORKSPACE_MB=${CPPAW_CUBLAS_FP64_WORKSPACE_MB:-8192}
export OZAKI_CASE_PREFIX=${OZAKI_CASE_PREFIX:-gpu_recommended_ozaki}
export GPU_CASES=${GPU_CASES:-"${OZAKI_CASE_PREFIX}_native ${OZAKI_CASE_PREFIX}_dgemm ${OZAKI_CASE_PREFIX}_zgemm ${OZAKI_CASE_PREFIX}_zherk ${OZAKI_CASE_PREFIX}_all"}
export ONE_RANK_CPU_CASES=${ONE_RANK_CPU_CASES-"cpu nvhpc_cpu"}
export CPU_CASES=${CPU_CASES-"cpu nvhpc_cpu"}

"${HERE}/run_followup.sh"

validation_log="${FOLLOWUP_ROOT}/ozaki_validation.log"
: > "${validation_log}"
for empty_bands in ${EMPTY_BANDS_LIST}; do
  suite="${FOLLOWUP_ROOT}/empty${empty_bands}_nstep1_${GPU_RANKS}rank_gpu"
  reference="${suite}/${OZAKI_CASE_PREFIX}_native/rep01"
  for kernel in dgemm zgemm zherk; do
    candidate="${suite}/${OZAKI_CASE_PREFIX}_${kernel}/rep01"
    python3 "${HERE}/ozaki_validate.py" "${reference}" "${candidate}" \
      --expected-kpoints 1 | tee -a "${validation_log}"
  done
  python3 "${HERE}/ozaki_validate.py" "${reference}" \
    "${suite}/${OZAKI_CASE_PREFIX}_all/rep01" \
    --expected-kpoints 1 | tee -a "${validation_log}"
done

correctness_root="${FOLLOWUP_ROOT}/si2_force_stress_kpoints"
TEST=si2_ozaki \
CNTL_FILE="${HERE}/si2_ozaki.cntl" \
STRC_FILE="${HERE}/si2_ozaki.strc" \
NSTEPS=1 RANKS=1 REPEATS=1 EXPECTED_ENERGY=none \
CASES="${OZAKI_CASE_PREFIX}_native ${OZAKI_CASE_PREFIX}_all" \
CPPAW_CUBLAS_FP64_STRATEGY=eager \
RUN_ROOT="${correctness_root}" \
"${HERE}/run_benchmark.sh"

python3 "${HERE}/ozaki_validate.py" \
  "${correctness_root}/${OZAKI_CASE_PREFIX}_native/rep01" \
  "${correctness_root}/${OZAKI_CASE_PREFIX}_all/rep01" \
  --expected-kpoints 2 --require-force --require-stress \
  --require-ozaki ZGEMM | tee -a "${validation_log}"

echo "Ozaki benchmark data: ${FOLLOWUP_ROOT}"
