#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

THISDIR=$(cd "$(dirname "$0")/../.." && pwd)
DEVICE=auto
JOBS=${CPPAW_INSTALL_JOBS:-16}
PREFIX=
TORCH_PREFIX=${CPPAW_TORCH_PREFIX:-}
DOWNLOAD_MODEL=false
CUDA_COMPILER=
CUDA_ARCH=
CUDA_ROOT=
TORCH_CUDA_ARCH=

usage() {
  cat <<USAGE
usage: $0 [--device auto|cpu|cuda] [--prefix DIR] [--torch-prefix DIR]
          [--jobs N] [--download-model]
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --device) DEVICE=$2; shift 2 ;;
    --prefix) PREFIX=$2; shift 2 ;;
    --torch-prefix) TORCH_PREFIX=$2; shift 2 ;;
    --jobs) JOBS=$2; shift 2 ;;
    --download-model) DOWNLOAD_MODEL=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ${DEVICE} = auto ]]; then
  if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    DEVICE=cuda
  else
    DEVICE=cpu
  fi
fi
case "${DEVICE}" in cpu|cuda) ;; *) usage >&2; exit 2 ;; esac

PYTHON=${CPPAW_SKALA_PYTHON:-$(command -v python3 || true)}
if [[ -z ${PYTHON} ]]; then
  echo "python3 is required to locate a LibTorch installation" >&2
  exit 1
fi
if [[ -z ${TORCH_PREFIX} ]]; then
  TORCH_PREFIX=$(${PYTHON} -c 'import torch; print(torch.utils.cmake_prefix_path)' 2>/dev/null || true)
fi
if [[ -z ${TORCH_PREFIX} ]]; then
  echo "PyTorch/LibTorch was not found." >&2
  echo "Install PyTorch or pass --torch-prefix /path/to/torch." >&2
  exit 1
fi
if [[ ${DEVICE} = cuda ]]; then
  TORCH_CUDA=$(${PYTHON} -c 'import torch; print(int(torch.cuda.is_available()))' 2>/dev/null || echo 0)
  if [[ ${TORCH_CUDA} != 1 ]]; then
    echo "CUDA was detected, but the selected PyTorch installation has no usable CUDA support." >&2
    exit 1
  fi

  TORCH_CUDA_VERSION=$(${PYTHON} -c 'import torch; print(torch.version.cuda or "")')
  if [[ -z ${TORCH_CUDA_VERSION} ]]; then
    echo "The selected PyTorch installation does not report a CUDA toolkit version." >&2
    exit 1
  fi

  nvcc_version() {
    "$1" --version 2>/dev/null \
      | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' \
      | head -1
  }

  if [[ -n ${CUDACXX:-} ]]; then
    CUDA_COMPILER=$(command -v "${CUDACXX}" || true)
    if [[ -z ${CUDA_COMPILER} ]]; then
      echo "CUDACXX=${CUDACXX} is not executable." >&2
      exit 1
    fi
    if [[ $(nvcc_version "${CUDA_COMPILER}") != "${TORCH_CUDA_VERSION}" ]]; then
      echo "CUDACXX=${CUDA_COMPILER} does not match PyTorch CUDA ${TORCH_CUDA_VERSION}." >&2
      exit 1
    fi
  else
    NVCC_CANDIDATES=()
    NVCC_ON_PATH=$(command -v nvcc || true)
    [[ -n ${NVCC_ON_PATH} ]] && NVCC_CANDIDATES+=("${NVCC_ON_PATH}")
    [[ -x /usr/local/cuda-${TORCH_CUDA_VERSION}/bin/nvcc ]] \
      && NVCC_CANDIDATES+=("/usr/local/cuda-${TORCH_CUDA_VERSION}/bin/nvcc")
    [[ -x /usr/local/cuda/bin/nvcc ]] \
      && NVCC_CANDIDATES+=("/usr/local/cuda/bin/nvcc")
    if [[ -n ${NVHPC_ROOT:-} \
          && -x ${NVHPC_ROOT}/cuda/${TORCH_CUDA_VERSION}/bin/nvcc ]]; then
      NVCC_CANDIDATES+=("${NVHPC_ROOT}/cuda/${TORCH_CUDA_VERSION}/bin/nvcc")
    fi
    shopt -s nullglob
    for NVCC in /opt/nvidia/hpc_sdk/*/*/cuda/${TORCH_CUDA_VERSION}/bin/nvcc; do
      NVCC_CANDIDATES+=("${NVCC}")
    done
    shopt -u nullglob
    for NVCC in "${NVCC_CANDIDATES[@]}"; do
      if [[ -x ${NVCC} && $(nvcc_version "${NVCC}") = "${TORCH_CUDA_VERSION}" ]]; then
        CUDA_COMPILER=${NVCC}
        break
      fi
    done
    if [[ -z ${CUDA_COMPILER} ]]; then
      echo "No nvcc matching PyTorch CUDA ${TORCH_CUDA_VERSION} was found." >&2
      echo "Install the matching CUDA toolkit or set CUDACXX explicitly." >&2
      exit 1
    fi
  fi

  CUDA_ARCH=${CPPAW_SKALA_CUDA_ARCH:-$(${PYTHON} -c '
import torch
arches = torch.cuda.get_arch_list()
if arches:
    print(arches[0].removeprefix("sm_"))
else:
    major, minor = torch.cuda.get_device_capability()
    print(f"{major}{minor}")
')}
  CUDA_ARCH=${CUDA_ARCH#sm_}
  if [[ ${CUDA_ARCH} =~ ^([0-9]+)([0-9])([a-z]?)$ ]]; then
    TORCH_CUDA_ARCH="${BASH_REMATCH[1]}.${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
  elif [[ ${CUDA_ARCH} =~ ^([0-9]+)\.([0-9])([a-z]?)$ ]]; then
    TORCH_CUDA_ARCH="${BASH_REMATCH[1]}.${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
    CUDA_ARCH="${BASH_REMATCH[1]}${BASH_REMATCH[2]}${BASH_REMATCH[3]}"
  else
    echo "Invalid CUDA architecture '${CUDA_ARCH}'. Use a value such as 80 or 8.0." >&2
    exit 1
  fi
  CUDA_ROOT=$(cd "$(dirname "${CUDA_COMPILER}")/.." && pwd)
fi

if [[ -z ${PREFIX} ]]; then
  PREFIX=${THISDIR}/bin/skala_ftorch_${DEVICE}
fi
BUILD=${THISDIR}/bin/Build_skala_ftorch_${DEVICE}
FC=${FC:-$(command -v nvfortran || command -v gfortran || true)}
CXX=${CXX:-$(command -v g++ || command -v clang++ || true)}
DEVICE_CMAKE=$(printf '%s' "${DEVICE}" | tr '[:lower:]' '[:upper:]')
if [[ -z ${FC} || -z ${CXX} ]]; then
  echo "a Fortran compiler and a C++17 compiler are required" >&2
  exit 1
fi

CUDA_CMAKE_ARGS=()
if [[ ${DEVICE} = cuda ]]; then
  CUDA_CMAKE_ARGS+=("-DCMAKE_CUDA_COMPILER=${CUDA_COMPILER}")
  CUDA_CMAKE_ARGS+=("-DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCH}")
  CUDA_CMAKE_ARGS+=("-DCUDAToolkit_ROOT=${CUDA_ROOT}")
  CUDA_CMAKE_ARGS+=("-DCUDA_TOOLKIT_ROOT_DIR=${CUDA_ROOT}")
  CUDA_CMAKE_ARGS+=("-DCUDA_NVCC_EXECUTABLE=${CUDA_COMPILER}")
  export CUDA_HOME=${CUDA_ROOT}
  export CUDA_PATH=${CUDA_ROOT}
  export PATH="${CUDA_ROOT}/bin:${PATH}"
  export TORCH_CUDA_ARCH_LIST=${TORCH_CUDA_ARCH}
  echo "Using nvcc ${CUDA_COMPILER} for PyTorch CUDA ${TORCH_CUDA_VERSION}"
  echo "Using CMake CUDA architecture ${CUDA_ARCH}"
fi

cmake -S "${THISDIR}/src/SkalaBridge" -B "${BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${TORCH_PREFIX}" \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCPPAW_SKALA_DEVICE="${DEVICE_CMAKE}" \
  "${CUDA_CMAKE_ARGS[@]}"
cmake --build "${BUILD}" --parallel "${JOBS}"
cmake --install "${BUILD}"

MODEL=
if [[ ${DOWNLOAD_MODEL} = true ]]; then
  MODEL=$(${THISDIR}/src/Buildtools/paw_skala_model.sh "${DEVICE}" "${PREFIX}/models")
  "${PREFIX}/bin/cppaw_skala_smoke" "${MODEL}" "${DEVICE}"
fi

echo "CP-PAW Skala bridge installed in ${PREFIX}"
if [[ -n ${MODEL} ]]; then
  echo "Skala model and bridge smoke test passed with ${MODEL}"
fi
