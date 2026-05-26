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
NSTEPS=${NSTEPS:-1}
RANKS=${RANKS:-1}
TIMEOUT=${TIMEOUT:-600}
RUN_ROOT=${RUN_ROOT:-"${HERE}/runs/${TEST}-nsys-nstep${NSTEPS}-${RANKS}ranks-$(date +%Y%m%d-%H%M%S)"}
MPI_ARGS=${MPI_ARGS:---mca coll ^hcoll}
EXE=${EXE:-"${ROOT}/bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x"}
NSYS_TRACE=${NSYS_TRACE:-cuda,nvtx,osrt}
NSYS_STATS=${NSYS_STATS:-true}
NSYS_MPI_MODE=${NSYS_MPI_MODE:-outer}

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

if ! command -v nsys >/dev/null 2>&1; then
  echo "nsys was not found in PATH." >&2
  exit 1
fi
if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${EXE}" >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}"
if [[ ! -f "${CNTL_FILE}" ]]; then
  echo "Control file not found: ${CNTL_FILE}" >&2
  exit 1
fi
if [[ ! -f "${STRC_FILE}" ]]; then
  echo "Structure file not found: ${STRC_FILE}" >&2
  exit 1
fi
cp "${CNTL_FILE}" "${RUN_ROOT}/${TEST}.cntl"
cp "${STRC_FILE}" "${RUN_ROOT}/${TEST}.strc"
cp "${HERE}/profile_summary.py" "${RUN_ROOT}/"
cp "${ROOT}/tests/fulltests/si2/stp.cntl" "${RUN_ROOT}/"
perl -0pi -e "s/NSTEP\\s*=\\s*\\d+/NSTEP=${NSTEPS}/" "${RUN_ROOT}/${TEST}.cntl"
if [[ -n "${EMPTY_BANDS:-}" ]]; then
  perl -0pi -e "s/EMPTY\\s*=\\s*\\d+/EMPTY=${EMPTY_BANDS}/" "${RUN_ROOT}/${TEST}.strc"
fi

cd "${RUN_ROOT}"
echo "${RUN_ROOT}" > "${HERE}/runs/latest_nsys"

if [[ "${RANKS}" -gt 1 && "${NSYS_MPI_MODE}" == "per_rank" ]]; then
  cat > nsys_rank_wrapper.sh <<EOF
#!/usr/bin/env bash
set -euo pipefail
rank=\${OMPI_COMM_WORLD_RANK:-\${PMI_RANK:-\${SLURM_PROCID:-0}}}
export CPPAW_ACCEL_PROFILE_FILE=nsys_profile
exec nsys profile --force-overwrite=true --stats="${NSYS_STATS}" --trace="${NSYS_TRACE}" \\
  -o "nsys_rank\${rank}" "${EXE}" "${TEST}.cntl"
EOF
  chmod +x nsys_rank_wrapper.sh
  # shellcheck disable=SC2086
  CMD="${MPIRUN} ${MPI_ARGS} -np ${RANKS} ./nsys_rank_wrapper.sh"
elif [[ "${RANKS}" -gt 1 ]]; then
  # CP-PAW exits MPI runs via MPI_ABORT(0), so keep nsys outside the MPI job
  # by default to let it flush the combined report.
  CMD="env CPPAW_ACCEL_PROFILE_FILE=nsys_profile nsys profile --force-overwrite=true --stats=${NSYS_STATS} --trace=${NSYS_TRACE} -o nsys_mpi ${MPIRUN} ${MPI_ARGS} -np ${RANKS} ${EXE} ${TEST}.cntl"
else
  CMD="env CPPAW_ACCEL_PROFILE_FILE=nsys_profile nsys profile --force-overwrite=true --stats=${NSYS_STATS} --trace=${NSYS_TRACE} -o nsys ${EXE} ${TEST}.cntl"
fi

echo "running Nsight Systems: ${CMD}"
# shellcheck disable=SC2086
/usr/bin/time -p timeout "${TIMEOUT}s" ${CMD} > out.log 2> err.log

python3 profile_summary.py nsys_profile*.csv > summary.txt
grep -E "^real |^user |^sys " err.log > time.txt || true
echo "Nsight data: ${RUN_ROOT}"
