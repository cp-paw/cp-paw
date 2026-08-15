#!/bin/sh
set -eu

atom=${1:-2}
axis=${2:-1}
case "$atom" in *[!0-9]*|'') echo "usage: $0 [atom [axis]]" >&2; exit 2 ;; esac
case "$axis" in 1|2|3) ;; *) echo "usage: $0 [atom [axis]]" >&2; exit 2 ;; esac
test "$atom" -gt 0 || { echo "atom must be positive" >&2; exit 2; }

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
step=${SKALA_FORCE_FD_STEP:-0.0003}
tolerance=${SKALA_FORCE_FD_TOLERANCE:-0.002}
ranks=${MPI_RANKS:-1}
mpiexec=${MPIEXEC:-mpirun}
work=$(mktemp -d "${TMPDIR:-/tmp}/cppaw-skala-force.XXXXXX")
keep=${SKALA_FORCE_FD_KEEP:-0}
cleanup() {
  if test "$keep" = 1; then
    echo "force finite-difference files retained in $work"
  else
    rm -rf "$work"
  fi
}
trap cleanup EXIT HUP INT TERM

cp "$here/../si2/stp.cntl" "$work/stp.cntl"
ln -s "$SKALA_MODEL" "$work/model.fun"

minus_step=$(awk -v h="$step" 'BEGIN { printf "%.17g", -h }')

make_control() {
  name=$1
  sed -e "s|NAME='../si2/stp.cntl'|NAME='stp.cntl'|" \
      -e 's/START=T/START=F/' \
      -e 's/!CELL MOVE=T/!CELL MOVE=F/' \
      "$here/skala_si2.cntl" >"$work/$name.cntl"
  cp "$SKALA_STRUCTURE" "$work/$name.strc"
}

run_case() {
  name=$1
  delta=$2
  make_control "$name"
  "$python" "$here/displace_restart.py" \
    "$SKALA_RESTART" "$work/$name.rstrt" "$atom" "$axis" "$delta"
  if test "$ranks" -eq 1; then
    (cd "$work" && "$PAWX" "$name.cntl" >"$name.out" 2>"$name.err")
  else
    (cd "$work" && "$mpiexec" -np "$ranks" "$PAWX" "$name.cntl" \
      >"$name.out" 2>"$name.err")
  fi
  grep -q "PROGRAM FINISHED" "$work/$name.prot"
}

run_case base 0
run_case plus "$step"
run_case minus "$minus_step"

energy() {
  awk '
    /SKALA TOTAL FORCE DIAGNOSTIC/ { section = 1; next }
    section && /^TOTAL ENERGY/ { value = $NF }
    END { if (value == "") exit 1; print value }
  ' "$1"
}

eplus=$(energy "$work/plus.prot")
eminus=$(energy "$work/minus.prot")
finite_difference=$(awk -v ep="$eplus" -v em="$eminus" -v h="$step" \
  'BEGIN { printf "%.14e", -(ep - em) / (2.0 * h) }')

analytic=$(awk -v atom="$atom" -v axis="$axis" '
  /SKALA TOTAL FORCE DIAGNOSTIC/ { section = 1; next }
  section && /^ATOM[[:space:]]/ && $2 == atom { value = $(2 + axis) }
  END { if (value == "") exit 1; printf "%.14e", value }
' "$work/base.prot")

awk -v atom="$atom" -v axis="$axis" -v h="$step" \
    -v analytic="$analytic" -v fd="$finite_difference" \
    -v tolerance="$tolerance" '
  BEGIN {
    difference = fd - analytic
    absolute = difference < 0.0 ? -difference : difference
    printf "SKALA FORCE FD atom=%d axis=%d step=%s analytic=% .14e fd=% .14e difference=% .6e\n", \
           atom, axis, h, analytic, fd, difference
    if (absolute > tolerance) {
      printf "TEST FAILED: absolute force error %.6e exceeds %.6e\n", \
             absolute, tolerance > "/dev/stderr"
      exit 1
    }
  }
'
echo "TEST PASSED"
