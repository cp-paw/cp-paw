#!/bin/bash
set -euo pipefail
output=${1:?Usage: run_isolated.sh new-output-directory}
: "${PAWX:?PAWX must name a Skala-enabled executable}"
: "${SKALA_MODEL:?SKALA_MODEL must name a Skala 1.1 model}"
device=${SKALA_DEVICE:-AUTO}
radial=${SKALA_RADIAL_POINTS:-200}
angular=${SKALA_LEBEDEV_EXACTNESS:-53}
case "$device" in CPU|CUDA|AUTO) ;; *) echo 'Invalid SKALA_DEVICE' >&2; exit 2 ;; esac
for value in "$radial" "$angular"; do
  case "$value" in *[!0-9]*|'') echo 'Invalid quadrature size' >&2; exit 2 ;; esac
  test "$value" -gt 0 || exit 2
done
PAWX=$(realpath "$PAWX")
SKALA_MODEL=$(realpath "$SKALA_MODEL")
here=$(cd "$(dirname "$0")" && pwd)
# CP-PAW appends to existing protocol files; each test needs a fresh directory.
mkdir "$output"
mkdir "$output/si2" "$output/skala_si2"
cp "$here/../si2/stp.cntl" "$output/si2/"
cp "$here/skala_si2.strc" "$here/analyse.sh" "$here/Makefile" "$output/skala_si2/"
sed -e "s/DEVICE='AUTO'/DEVICE='$device'/" \
    -e "s/RADIALPOINTS=200/RADIALPOINTS=$radial/" \
    -e "s/LEBEDEVEXACTNESS=53/LEBEDEVEXACTNESS=$angular/" \
    "$here/skala_si2.cntl" > "$output/skala_si2/skala_si2.cntl"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
"${MAKE:-make}" -C "$output/skala_si2" all PAWX="$PAWX" SKALA_MODEL="$SKALA_MODEL"
