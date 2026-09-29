#!/bin/bash
set -euo pipefail
output=${1:?Usage: run_isolated.sh new-output-directory}
: "${PAWX:?PAWX must name a Skala-enabled executable}"
: "${SKALA_MODEL:?SKALA_MODEL must name a Skala 1.1 model}"
device=${SKALA_DEVICE:-AUTO}
radial=${SKALA_RADIAL_POINTS:-200}
angular=${SKALA_LEBEDEV_EXACTNESS:-53}
orientations=${SKALA_LEBEDEV_ORIENTATIONS:-1}
shells=${SKALA_IMAGE_SHELLS:-1}
cutoff=${SKALA_CUTOFF:-20}
cdual=${SKALA_CDUAL:-2}
case "$device" in CPU|CUDA|AUTO) ;; *) echo 'Invalid SKALA_DEVICE' >&2; exit 2 ;; esac
for value in "$radial" "$angular" "$orientations" "$shells"; do
  case "$value" in *[!0-9]*|'') echo 'Invalid quadrature size' >&2; exit 2 ;; esac
  test "$value" -gt 0 || exit 2
done
test "$angular" -le 53 || { echo 'Maximum supported Lebedev exactness is 53' >&2; exit 2; }
for value in "$cutoff" "$cdual"; do
  case "$value" in *[!0-9.eE+-]*|'') echo 'Invalid Fourier cutoff' >&2; exit 2 ;; esac
  awk -v value="$value" 'BEGIN { exit !(value + 0 > 0) }' || exit 2
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
    -e "s/LEBEDEVORIENTATIONS=1/LEBEDEVORIENTATIONS=$orientations IMAGESHELLS=$shells/" \
    -e "s/EPWPSI=20 CDUAL=2/EPWPSI=$cutoff CDUAL=$cdual/" \
    "$here/skala_si2.cntl" > "$output/skala_si2/skala_si2.cntl"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
"${MAKE:-make}" -C "$output/skala_si2" all PAWX="$PAWX" SKALA_MODEL="$SKALA_MODEL"
