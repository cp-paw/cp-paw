#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

DEVICE=${1:-cpu}
DESTINATION=${2:-"${PWD}/models"}

case "${DEVICE}" in
  cpu)
    FILE=skala-1.1-rev1.fun
    SHA256=7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd
    ;;
  cuda)
    FILE=skala-1.1-rev1-cuda.fun
    SHA256=f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223
    ;;
  *)
    echo "usage: $0 cpu|cuda [destination]" >&2
    exit 2
    ;;
esac

mkdir -p "${DESTINATION}"
TARGET=${DESTINATION}/${FILE}
URL=https://huggingface.co/microsoft/skala-1.1/resolve/main/${FILE}

if [[ ! -f ${TARGET} ]]; then
  rm -f "${TARGET}.tmp"
  DOWNLOADED=false
  if command -v curl >/dev/null 2>&1 \
      && curl --fail --location --retry 3 --output "${TARGET}.tmp" "${URL}"; then
    DOWNLOADED=true
  elif command -v wget >/dev/null 2>&1 \
      && wget --output-document="${TARGET}.tmp" "${URL}"; then
    DOWNLOADED=true
  elif command -v python3 >/dev/null 2>&1 \
      && python3 - "${URL}" "${TARGET}.tmp" <<'PY'
import sys
import urllib.request

urllib.request.urlretrieve(sys.argv[1], sys.argv[2])
PY
  then
    DOWNLOADED=true
  fi
  if [[ ${DOWNLOADED} != true ]]; then
    echo "could not download the Skala model with curl, wget, or python3" >&2
    exit 1
  fi
  mv "${TARGET}.tmp" "${TARGET}"
fi

if command -v sha256sum >/dev/null 2>&1; then
  ACTUAL=$(sha256sum "${TARGET}" | awk '{print $1}')
else
  ACTUAL=$(shasum -a 256 "${TARGET}" | awk '{print $1}')
fi
if [[ ${ACTUAL} != ${SHA256} ]]; then
  echo "Skala model hash mismatch for ${TARGET}" >&2
  echo "expected ${SHA256}" >&2
  echo "actual   ${ACTUAL}" >&2
  exit 1
fi

printf '%s\n' "${TARGET}"
