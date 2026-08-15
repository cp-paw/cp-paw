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
radial_list=${SKALA_RADIAL_POINT_LIST:-"100 200 400"}
lebedev_list=${SKALA_LEBEDEV_EXACTNESS_LIST:-"17"}
orientation_list=${SKALA_LEBEDEV_ORIENTATION_LIST:-"1"}
device=${SKALA_DEVICE:-AUTO}
ranks=${MPI_RANKS:-1}
mpiexec=${MPIEXEC:-mpirun}
case "$device" in AUTO|CPU|CUDA) ;; *) echo "SKALA_DEVICE must be AUTO, CPU, or CUDA" >&2; exit 2 ;; esac
work=$(mktemp -d "${TMPDIR:-/tmp}/cppaw-skala-quadrature.XXXXXX")
keep=${SKALA_QUADRATURE_KEEP:-0}
cleanup() {
  if test "$keep" = 1; then
    echo "quadrature-convergence files retained in $work"
  else
    find "$work" -depth -delete
  fi
}
trap cleanup EXIT HUP INT TERM

case "$ranks" in *[!0-9]*|'') echo "MPI_RANKS must be an integer" >&2; exit 2 ;; esac
test "$ranks" -gt 0 || { echo "MPI_RANKS must be positive" >&2; exit 2; }
for value in $radial_list; do
  case "$value" in *[!0-9]*|'') echo "SKALA_RADIAL_POINT_LIST must contain integers" >&2; exit 2 ;; esac
  test "$value" -gt 0 || { echo "radial point counts must be positive" >&2; exit 2; }
done
for value in $lebedev_list; do
  case "$value" in *[!0-9]*|'') echo "SKALA_LEBEDEV_EXACTNESS_LIST must contain integers" >&2; exit 2 ;; esac
  test "$value" -gt 0 || { echo "Lebedev exactness values must be positive" >&2; exit 2; }
done
for value in $orientation_list; do
  case "$value" in *[!0-9]*|'') echo "SKALA_LEBEDEV_ORIENTATION_LIST must contain integers" >&2; exit 2 ;; esac
  test "$value" -gt 0 || { echo "Lebedev orientation counts must be positive" >&2; exit 2; }
done

cp "$here/../si2/stp.cntl" "$work/stp.cntl"
ln -s "$SKALA_MODEL" "$work/model.fun"
printf "radial\tlebedev\torientations\trows\telectrons\tmodel_xc\n" >"$work/raw.tsv"

run_case() {
  radial=$1
  lebedev=$2
  orientations=$3
  name="radial${radial}_lebedev${lebedev}_orient${orientations}"
  sed -e "s|NAME='../si2/stp.cntl'|NAME='stp.cntl'|" \
      -e 's/START=T/START=F/' \
      -e 's/!CELL MOVE=T/!CELL MOVE=F/' \
      -e 's/CHECK=T/CHECK=F/' \
      -e "s/RADIALPOINTS=200/RADIALPOINTS=$radial/" \
      -e "s/LEBEDEVEXACTNESS=53/LEBEDEVEXACTNESS=$lebedev/" \
      -e "s/LEBEDEVORIENTATIONS=1/LEBEDEVORIENTATIONS=$orientations/" \
      -e "s/DEVICE='AUTO'/DEVICE='$device'/" \
      "$here/skala_si2.cntl" >"$work/$name.cntl"
  cp "$SKALA_STRUCTURE" "$work/$name.strc"
  cp "$SKALA_RESTART" "$work/$name.rstrt"
  if test "$ranks" -eq 1; then
    (cd "$work" && "$PAWX" "$name.cntl" >"$name.out" 2>"$name.err")
  else
    (cd "$work" && "$mpiexec" -np "$ranks" "$PAWX" "$name.cntl" \
      >"$name.out" 2>"$name.err")
  fi
  grep -q "PROGRAM FINISHED" "$work/$name.prot"
  awk -v radial="$radial" -v lebedev="$lebedev" \
      -v orientations="$orientations" '
    /^HYBRID-GRID ROWS/ { rows = $NF }
    /^COMPOSITE ELECTRONS/ { electrons = $NF }
    /^MODEL XC ENERGY/ { model_xc = $NF }
    END {
      if (rows == "" || electrons == "" || model_xc == "") exit 1
      printf "%d\t%d\t%d\t%d\t%.14e\t%.14e\n", radial, lebedev, orientations, rows, electrons, model_xc
    }
  ' "$work/$name.prot" >>"$work/raw.tsv"
}

for orientations in $orientation_list; do
  for lebedev in $lebedev_list; do
    for radial in $radial_list; do
      run_case "$radial" "$lebedev" "$orientations"
    done
  done
done

awk -F '\t' '
  NR == 1 { next }
  {
    count++
    radial[count] = $1
    lebedev[count] = $2
    orientations[count] = $3
    rows[count] = $4
    electrons[count] = $5
    model_xc[count] = $6
  }
  END {
    if (count == 0) exit 1
    print "RADIAL LEBEDEV ORIENTATIONS ROWS ELECTRONS DELTA_ELECTRONS MODEL_XC DELTA_MODEL_XC"
    for (i = 1; i <= count; i++)
      printf "%6d %7d %12d %8d % .12e % .6e % .12e % .6e\n", \
             radial[i], lebedev[i], orientations[i], rows[i], electrons[i], \
             electrons[i] - electrons[count], model_xc[i], \
             model_xc[i] - model_xc[count]
  }
' "$work/raw.tsv"

echo "TEST PASSED"
