#!/bin/sh
set -eu

case_name=${1:-isotropic}
case "$case_name" in
  isotropic|xx|xy) ;;
  *) echo "usage: $0 [isotropic|xx|xy]" >&2; exit 2 ;;
esac

: "${PAWX:?PAWX must name a CP-PAW executable with Skala support}"
: "${SKALA_MODEL:?SKALA_MODEL must name a Skala TorchScript model}"
: "${SKALA_RESTART:?SKALA_RESTART must name a converged CP-PAW restart}"
: "${SKALA_STRUCTURE:?SKALA_STRUCTURE must match SKALA_RESTART}"

for file in "$PAWX" "$SKALA_MODEL" "$SKALA_RESTART" "$SKALA_STRUCTURE"; do
  test -e "$file" || { echo "missing input: $file" >&2; exit 2; }
done
PAWX=$(realpath "$PAWX")
SKALA_MODEL=$(realpath "$SKALA_MODEL")
SKALA_RESTART=$(realpath "$SKALA_RESTART")
SKALA_STRUCTURE=$(realpath "$SKALA_STRUCTURE")

here=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
python=${PYTHON:-python3}
step=${SKALA_STRESS_FD_STEP:-0.0003}
tolerance=${SKALA_STRESS_FD_TOLERANCE:-0.002}
radial_points=${SKALA_RADIAL_POINTS:-200}
lebedev_exactness=${SKALA_LEBEDEV_EXACTNESS:-53}
lebedev_orientations=${SKALA_LEBEDEV_ORIENTATIONS:-1}
device=${SKALA_DEVICE:-AUTO}
ranks=${MPI_RANKS:-1}
mpiexec=${MPIEXEC:-mpirun}
case "$device" in AUTO|CPU|CUDA) ;; *) echo "SKALA_DEVICE must be AUTO, CPU, or CUDA" >&2; exit 2 ;; esac
case "$radial_points" in *[!0-9]*|'') echo "SKALA_RADIAL_POINTS must be an integer" >&2; exit 2 ;; esac
case "$lebedev_exactness" in *[!0-9]*|'') echo "SKALA_LEBEDEV_EXACTNESS must be an integer" >&2; exit 2 ;; esac
case "$lebedev_orientations" in *[!0-9]*|'') echo "SKALA_LEBEDEV_ORIENTATIONS must be an integer" >&2; exit 2 ;; esac
test "$radial_points" -gt 0 || { echo "SKALA_RADIAL_POINTS must be positive" >&2; exit 2; }
test "$lebedev_exactness" -gt 0 || { echo "SKALA_LEBEDEV_EXACTNESS must be positive" >&2; exit 2; }
test "$lebedev_orientations" -gt 0 || { echo "SKALA_LEBEDEV_ORIENTATIONS must be positive" >&2; exit 2; }
work=$(mktemp -d "${TMPDIR:-/tmp}/cppaw-skala-stress.XXXXXX")
keep=${SKALA_STRESS_FD_KEEP:-0}
cleanup() {
  if test "$keep" = 1; then
    echo "stress finite-difference files retained in $work"
  else
    rm -rf "$work"
  fi
}
trap cleanup EXIT HUP INT TERM

cp "$here/../si2/stp.cntl" "$work/stp.cntl"
ln -s "$SKALA_MODEL" "$work/model.fun"

half_step=$(awk -v h="$step" 'BEGIN { printf "%.17g", 0.5 * h }')
minus_step=$(awk -v h="$step" 'BEGIN { printf "%.17g", -h }')
minus_half=$(awk -v h="$step" 'BEGIN { printf "%.17g", -0.5 * h }')

case "$case_name" in
  isotropic)
    plus="$step 0 0 0 $step 0 0 0 $step"
    minus="$minus_step 0 0 0 $minus_step 0 0 0 $minus_step"
    ;;
  xx)
    plus="$step 0 0 0 0 0 0 0 0"
    minus="$minus_step 0 0 0 0 0 0 0 0"
    ;;
  xy)
    plus="0 $half_step 0 $half_step 0 0 0 0 0"
    minus="0 $minus_half 0 $minus_half 0 0 0 0 0"
    ;;
esac

make_control() {
  name=$1
  sed -e "s|NAME='../si2/stp.cntl'|NAME='stp.cntl'|" \
      -e 's/START=T/START=F/' \
      -e "s/RADIALPOINTS=200/RADIALPOINTS=$radial_points/" \
      -e "s/LEBEDEVEXACTNESS=53/LEBEDEVEXACTNESS=$lebedev_exactness/" \
      -e "s/LEBEDEVORIENTATIONS=1/LEBEDEVORIENTATIONS=$lebedev_orientations/" \
      -e "s/DEVICE='AUTO'/DEVICE='$device'/" \
      "$here/skala_si2.cntl" >"$work/$name.cntl"
  cp "$SKALA_STRUCTURE" "$work/$name.strc"
}

run_case() {
  name=$1
  shift
  make_control "$name"
  "$python" "$here/deform_restart.py" \
    "$SKALA_RESTART" "$work/$name.rstrt" "$@"
  if test "$ranks" -eq 1; then
    (cd "$work" && "$PAWX" "$name.cntl" >"$name.out" 2>"$name.err")
  else
    (cd "$work" && "$mpiexec" -np "$ranks" "$PAWX" "$name.cntl" \
      >"$name.out" 2>"$name.err")
  fi
  grep -q "PROGRAM FINISHED" "$work/$name.prot"
}

run_case base 0 0 0 0 0 0 0 0 0
# shellcheck disable=SC2086
run_case plus $plus
# shellcheck disable=SC2086
run_case minus $minus

energy() {
  awk '
    /SKALA TOTAL FORCE DIAGNOSTIC/ { section = 1; next }
    section && /^TOTAL ENERGY[[:space:]]+[-+0-9]/ { value = $NF; exit }
    END { if (value == "") exit 1; print value }
  ' "$1"
}

eplus=$(energy "$work/plus.prot")
eminus=$(energy "$work/minus.prot")
finite_difference=$(awk -v ep="$eplus" -v em="$eminus" -v h="$step" \
  'BEGIN { printf "%.14e", (ep - em) / (2.0 * h) }')

analytic=$(awk -v case_name="$case_name" '
  /^TOTAL D E \/ D STRAIN/ {
    row++
    for (column = 1; column <= 3; column++) tensor[row, column] = $(NF - 3 + column)
  }
  END {
    if (row != 3) exit 1
    if (case_name == "isotropic") value = tensor[1,1] + tensor[2,2] + tensor[3,3]
    else if (case_name == "xx") value = tensor[1,1]
    else value = 0.5 * (tensor[1,2] + tensor[2,1])
    printf "%.14e", value
  }
' "$work/base.prot")

awk -v name="$case_name" -v h="$step" -v analytic="$analytic" \
    -v fd="$finite_difference" -v tolerance="$tolerance" '
  BEGIN {
    difference = fd - analytic
    absolute = difference < 0.0 ? -difference : difference
    printf "SKALA STRESS FD %-9s step=%s analytic=% .14e fd=% .14e difference=% .6e\n", \
           name, h, analytic, fd, difference
    if (absolute > tolerance) {
      printf "TEST FAILED: absolute stress error %.6e exceeds %.6e\n", \
             absolute, tolerance > "/dev/stderr"
      exit 1
    }
  }
'
echo "TEST PASSED"
