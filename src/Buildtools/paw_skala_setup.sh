#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

THISDIR=$(cd "$(dirname "$0")/../.." && pwd)
DEVICE=auto
JOBS=${CPPAW_INSTALL_JOBS:-16}
PREFIX=
TORCH_PREFIX=${CPPAW_TORCH_PREFIX:-}
DOWNLOAD_MODEL=false

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

cmake -S "${THISDIR}/src/SkalaBridge" -B "${BUILD}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${TORCH_PREFIX}" \
  -DCMAKE_Fortran_COMPILER="${FC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCPPAW_SKALA_DEVICE="${DEVICE_CMAKE}"
cmake --build "${BUILD}" --parallel "${JOBS}"
cmake --install "${BUILD}"

MODEL=
if [[ ${DOWNLOAD_MODEL} = true ]]; then
  MODEL=$(${THISDIR}/src/Buildtools/paw_skala_model.sh "${DEVICE}" "${PREFIX}/models")
fi

echo "CP-PAW Skala bridge installed in ${PREFIX}"
if [[ -n ${MODEL} ]]; then
  echo "Smoke test: ${PREFIX}/bin/cppaw_skala_smoke ${MODEL} ${DEVICE}"
fi
