#!/bin/bash
################################################################################
##  Run a command with NVIDIA NVBLAS interposition.
##
##  NVBLAS is a CUDA/cuBLASXt based drop-in BLAS layer. It intercepts selected
##  BLAS-3 routines such as GEMM, SYRK and HERK while falling back to a CPU BLAS
##  library for routines it does not handle.
################################################################################
set -euo pipefail

function usage {
  cat <<'EOF'
Usage:
  paw_nvblas.sh [--] command [args...]

Environment:
  CPPAW_NVBLAS_LIB             full path to libnvblas.so
  CPPAW_NVBLAS_CPU_BLAS_LIB    full path to the CPU BLAS fallback
  CPPAW_NVBLAS_GPU_LIST        NVBLAS GPU list, default ALL0
  CPPAW_NVBLAS_TILE_DIM        NVBLAS tile dimension, default 2048
  CPPAW_NVBLAS_AUTOPIN         yes/no, default no
  CPPAW_NVBLAS_TRACE           yes/no, default no
  CPPAW_NVBLAS_LOGFILE         NVBLAS log file, default nvblas.log
  NVBLAS_CONFIG_FILE           use an existing NVBLAS config instead of creating one
EOF
}

case "${1:-}" in
  -h|--help)
    usage
    exit 0
    ;;
  --)
    shift
    ;;
esac

if [[ $# -eq 0 ]] ; then
  usage >&2
  exit 2
fi

function platform_name {
  case "$(uname -s)_$(uname -m)" in
    Linux_aarch64|Linux_arm64) echo "Linux_aarch64" ;;
    Linux_x86_64) echo "Linux_x86_64" ;;
    *) echo "" ;;
  esac
}

function find_nvhpc_root {
  local platform
  local compiler_path
  local candidate

  if [[ -n "${NVHPC_ROOT:-}" && -x "${NVHPC_ROOT}/compilers/bin/nvfortran" ]] ; then
    echo "${NVHPC_ROOT}"
    return 0
  fi

  compiler_path=$(command -v nvfortran 2>/dev/null || true)
  if [[ -n ${compiler_path} ]] ; then
    dirname "$(dirname "${compiler_path}")"
    return 0
  fi

  platform=$(platform_name)
  if [[ -n ${platform} ]] ; then
    for candidate in /opt/nvidia/hpc_sdk/${platform}/* ; do
      if [[ -x ${candidate}/compilers/bin/nvfortran ]] ; then
        echo "${candidate}"
      fi
    done | sort -V | tail -1
  fi
}

function first_existing_file {
  local candidate
  for candidate in "$@" ; do
    if [[ -f ${candidate} ]] ; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

function find_nvblas {
  local root=$1
  local candidate

  if [[ -n "${CPPAW_NVBLAS_LIB:-}" ]] ; then
    first_existing_file "${CPPAW_NVBLAS_LIB}"
    return
  fi

  for candidate in \
      "${root}"/math_libs/*/targets/*/lib/libnvblas.so \
      "${root}"/math_libs/*/lib/libnvblas.so \
      /usr/local/cuda/lib64/libnvblas.so \
      /usr/local/cuda-*/lib64/libnvblas.so ; do
    if [[ -f ${candidate} ]] ; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

function find_cpu_blas {
  local root=$1
  local candidate

  if [[ -n "${CPPAW_NVBLAS_CPU_BLAS_LIB:-}" ]] ; then
    first_existing_file "${CPPAW_NVBLAS_CPU_BLAS_LIB}"
    return
  fi

  for candidate in \
      "${root}"/math_libs/nvpl/lib/libnvpl_blas_lp64_seq.so \
      "${root}"/math_libs/nvpl/lib/libnvpl_blas_lp64_gomp.so \
      /usr/lib*/libopenblas.so \
      /usr/lib*/libblas.so \
      /usr/local/lib*/libopenblas.so \
      /usr/local/lib*/libblas.so ; do
    if [[ -f ${candidate} ]] ; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

NVHPC_DETECTED_ROOT=$(find_nvhpc_root || true)
NVBLAS_LIB=$(find_nvblas "${NVHPC_DETECTED_ROOT}" || true)
CPU_BLAS_LIB=$(find_cpu_blas "${NVHPC_DETECTED_ROOT}" || true)

if [[ -z ${NVBLAS_LIB} ]] ; then
  echo "error in $0: libnvblas.so not found" >&2
  echo "set CPPAW_NVBLAS_LIB or NVHPC_ROOT" >&2
  exit 1
fi
if [[ -z ${CPU_BLAS_LIB} ]] ; then
  echo "error in $0: no CPU BLAS fallback found for NVBLAS" >&2
  echo "set CPPAW_NVBLAS_CPU_BLAS_LIB" >&2
  exit 1
fi

CREATED_CONFIG=
if [[ -z "${NVBLAS_CONFIG_FILE:-}" ]] ; then
  CREATED_CONFIG=$(mktemp "${TMPDIR:-/tmp}/cppaw-nvblas.XXXXXX.conf")
  export NVBLAS_CONFIG_FILE=${CREATED_CONFIG}
  {
    echo "NVBLAS_LOGFILE ${CPPAW_NVBLAS_LOGFILE:-nvblas.log}"
    if [[ "${CPPAW_NVBLAS_TRACE:-no}" = yes \
          || "${CPPAW_NVBLAS_TRACE:-no}" = true \
          || "${CPPAW_NVBLAS_TRACE:-no}" = 1 ]] ; then
      echo "NVBLAS_TRACE_LOG_ENABLED"
    fi
    echo "NVBLAS_CPU_BLAS_LIB ${CPU_BLAS_LIB}"
    echo "NVBLAS_GPU_LIST ${CPPAW_NVBLAS_GPU_LIST:-ALL0}"
    echo "NVBLAS_TILE_DIM ${CPPAW_NVBLAS_TILE_DIM:-2048}"
    if [[ "${CPPAW_NVBLAS_AUTOPIN:-no}" = yes \
          || "${CPPAW_NVBLAS_AUTOPIN:-no}" = true \
          || "${CPPAW_NVBLAS_AUTOPIN:-no}" = 1 ]] ; then
      echo "NVBLAS_AUTOPIN_MEM_ENABLED"
    fi
  } > "${NVBLAS_CONFIG_FILE}"
  chmod 600 "${NVBLAS_CONFIG_FILE}"
fi

NVBLAS_DIR=$(dirname "${NVBLAS_LIB}")
CPU_BLAS_DIR=$(dirname "${CPU_BLAS_LIB}")
export LD_LIBRARY_PATH="${NVBLAS_DIR}:${CPU_BLAS_DIR}:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="${NVBLAS_LIB}${LD_PRELOAD:+:${LD_PRELOAD}}"

"$@"
RC=$?

if [[ -n ${CREATED_CONFIG} ]] ; then
  rm -f "${CREATED_CONFIG}"
fi
exit ${RC}
