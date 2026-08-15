#!/bin/sh
set -eu

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
ranks=${MPI_RANKS:-4}
mpiexec=${MPIEXEC:-mpirun}
expected_kpoints=${SKALA_EXPECT_KPOINTS:-8}
energy_tolerance=${SKALA_MPI_ENERGY_TOLERANCE:-0.0000001}
force_tolerance=${SKALA_MPI_FORCE_TOLERANCE:-0.000002}
stress_tolerance=${SKALA_MPI_STRESS_TOLERANCE:-0.000002}
norm_tolerance=${SKALA_MPI_NORM_TOLERANCE:-0.00002}
radial_points=${SKALA_RADIAL_POINTS:-200}
lebedev_exactness=${SKALA_LEBEDEV_EXACTNESS:-53}
lebedev_orientations=${SKALA_LEBEDEV_ORIENTATIONS:-1}
work=$(mktemp -d "${TMPDIR:-/tmp}/cppaw-skala-mpi.XXXXXX")
keep=${SKALA_MPI_KEEP:-0}
case "$radial_points" in *[!0-9]*|'') echo "SKALA_RADIAL_POINTS must be an integer" >&2; exit 2 ;; esac
case "$lebedev_exactness" in *[!0-9]*|'') echo "SKALA_LEBEDEV_EXACTNESS must be an integer" >&2; exit 2 ;; esac
case "$lebedev_orientations" in *[!0-9]*|'') echo "SKALA_LEBEDEV_ORIENTATIONS must be an integer" >&2; exit 2 ;; esac
test "$radial_points" -gt 0 || { echo "SKALA_RADIAL_POINTS must be positive" >&2; exit 2; }
test "$lebedev_exactness" -gt 0 || { echo "SKALA_LEBEDEV_EXACTNESS must be positive" >&2; exit 2; }
test "$lebedev_orientations" -gt 0 || { echo "SKALA_LEBEDEV_ORIENTATIONS must be positive" >&2; exit 2; }
cleanup() {
  if test "$keep" = 1; then
    echo "MPI parity files retained in $work"
  else
    rm -rf "$work"
  fi
}
trap cleanup EXIT HUP INT TERM

case "$ranks" in *[!0-9]*|'') echo "MPI_RANKS must be an integer" >&2; exit 2 ;; esac
test "$ranks" -gt 1 || { echo "MPI_RANKS must be greater than one" >&2; exit 2; }

cp "$here/../si2/stp.cntl" "$work/stp.cntl"
ln -s "$SKALA_MODEL" "$work/model.fun"

make_case() {
  name=$1
  sed -e "s|NAME='../si2/stp.cntl'|NAME='stp.cntl'|" \
      -e 's/START=T/START=F/' \
      -e "s/RADIALPOINTS=200/RADIALPOINTS=$radial_points/" \
      -e "s/LEBEDEVEXACTNESS=53/LEBEDEVEXACTNESS=$lebedev_exactness/" \
      -e "s/LEBEDEVORIENTATIONS=1/LEBEDEVORIENTATIONS=$lebedev_orientations/" \
      "$here/skala_si2.cntl" >"$work/$name.cntl"
  cp "$SKALA_STRUCTURE" "$work/$name.strc"
  cp "$SKALA_RESTART" "$work/$name.rstrt"
}

make_case rank1
make_case rankn
(cd "$work" && "$mpiexec" -np 1 "$PAWX" rank1.cntl \
  >rank1.out 2>rank1.err)
(cd "$work" && "$mpiexec" -np "$ranks" "$PAWX" rankn.cntl \
  >rankn.out 2>rankn.err)
grep -q "PROGRAM FINISHED" "$work/rank1.prot"
grep -q "PROGRAM FINISHED" "$work/rankn.prot"

extract() {
  awk '
    /NUMBER OF K-POINTS/ { kpoints = $NF }
    /SKALA TOTAL FORCE DIAGNOSTIC/ { force_section = 1; next }
    force_section && /^TOTAL ENERGY[[:space:]]+[-+0-9]/ && energy == "" {
      energy = $NF
    }
    force_section && /^ATOM[[:space:]]/ {
      force_count++
      force_label[force_count] = "FORCE" $2
      for (i = 1; i <= 3; i++) force[force_count, i] = $(2 + i)
    }
    /^TOTAL D E \/ D STRAIN/ {
      stress_count++
      for (i = 1; i <= 3; i++) stress[stress_count, i] = $(NF - 3 + i)
    }
    /^SMOOTH SCALAR OPERATOR L2/ { scalar_norm = $NF }
    /^SMOOTH TAU OPERATOR L2/ { tau_norm = $NF }
    END {
      if (kpoints == "" || energy == "" || force_count == 0 || stress_count != 3 ||
          scalar_norm == "" || tau_norm == "") exit 1
      print "KPOINTS", kpoints
      print "ENERGY", energy
      print "SCALAR_NORM", scalar_norm
      print "TAU_NORM", tau_norm
      for (row = 1; row <= force_count; row++)
        print force_label[row], force[row, 1], force[row, 2], force[row, 3]
      for (row = 1; row <= stress_count; row++)
        print "STRESS" row, stress[row, 1], stress[row, 2], stress[row, 3]
    }
  ' "$1"
}

extract "$work/rank1.prot" >"$work/rank1.data"
extract "$work/rankn.prot" >"$work/rankn.data"

awk -v expected_kpoints="$expected_kpoints" \
    -v ranks="$ranks" \
    -v energy_tolerance="$energy_tolerance" \
    -v force_tolerance="$force_tolerance" \
    -v stress_tolerance="$stress_tolerance" \
    -v norm_tolerance="$norm_tolerance" '
  FNR == NR {
    reference[FNR] = $0
    reference_count = FNR
    next
  }
  {
    count = FNR
    left_count = split(reference[FNR], left)
    right_count = split($0, right)
    if (left_count != right_count || left[1] != right[1]) {
      print "TEST FAILED: MPI diagnostic layouts differ" > "/dev/stderr"
      failed = 1
      next
    }
    if (left[1] == "KPOINTS") {
      if (left[2] != expected_kpoints || right[2] != expected_kpoints) {
        print "TEST FAILED: unexpected k-point count" > "/dev/stderr"
        failed = 1
      }
      next
    }
    tolerance = energy_tolerance
    if (left[1] ~ /^FORCE/) tolerance = force_tolerance
    if (left[1] ~ /^STRESS/) tolerance = stress_tolerance
    if (left[1] ~ /_NORM$/) tolerance = norm_tolerance
    for (column = 2; column <= left_count; column++) {
      difference = right[column] - left[column]
      absolute = difference < 0.0 ? -difference : difference
      if (absolute > maximum[left[1]]) maximum[left[1]] = absolute
      if (absolute > tolerance) {
        printf "TEST FAILED: %s component %d difference %.6e exceeds %.6e\n", \
               left[1], column - 1, absolute, tolerance > "/dev/stderr"
        failed = 1
      }
    }
  }
  END {
    if (count != reference_count) {
      print "TEST FAILED: MPI diagnostic row counts differ" > "/dev/stderr"
      failed = 1
    }
    energy_max = maximum["ENERGY"] + 0.0
    force_max = 0.0
    stress_max = 0.0
    norm_max = 0.0
    for (label in maximum) {
      if (label ~ /^FORCE/ && maximum[label] > force_max) force_max = maximum[label]
      if (label ~ /^STRESS/ && maximum[label] > stress_max) stress_max = maximum[label]
      if (label ~ /_NORM$/ && maximum[label] > norm_max) norm_max = maximum[label]
    }
    printf "SKALA MPI PARITY ranks=1/%d energy=% .6e force=% .6e stress=% .6e norm=% .6e\n", \
           ranks, energy_max, force_max, stress_max, norm_max
    exit failed
  }
' "$work/rank1.data" "$work/rankn.data"
echo "TEST PASSED"
