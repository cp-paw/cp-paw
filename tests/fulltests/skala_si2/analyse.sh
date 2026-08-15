#!/bin/sh
set -eu

PROT=skala_si2.prot

value_for() {
  grep "$1" "$PROT" | tail -n 1 | awk '{print $NF}'
}

check_abs() {
  label=$1
  tolerance=$2
  value=$(value_for "$label")
  awk -v label="$label" -v value="$value" -v tolerance="$tolerance" '
    BEGIN {
      magnitude = value + 0.0
      if (magnitude < 0.0) magnitude = -magnitude
      if (value !~ /^[-+]?[0-9.]+[EeDd][-+]?[0-9]+$/ || magnitude >= tolerance) {
        printf "TEST FAILED: %s = %s (tolerance %s)\n", label, value, tolerance
        exit 1
      }
    }'
}

test "$(value_for 'NUMBER OF K-POINTS')" = "8" || {
  echo "TEST FAILED: expected eight k-points" >&2
  exit 1
}

check_abs "SKALA POSITIVE TAU CHECK" 1.0e-10
check_abs "GRAD ADJOINT DIFFERENCE" 1.0e-10
check_abs "TAU OPERATOR DIFFERENCE" 1.0e-10
check_abs "ONE-CENTER MATRIX DIFFERENCE" 1.0e-10

value_for "FULL ELECTRONIC OPERATOR CONTRACTION" >/dev/null
grep -q "PROGRAM FINISHED" "$PROT"
echo "TEST PASSED"
