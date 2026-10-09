#!/bin/sh
set -eu

PROT=skala_si2.prot

value_for() {
  grep "$1" "$PROT" | tail -n 1 | awk '{print $NF}'
}

check_abs() {
  label=$1
  tolerance=$2
  scale=${3:-1}
  awk -v label="$label" -v tolerance="$tolerance" -v scale="$scale" '
    index($0, label) == 1 && substr($0, length(label) + 1, 1) ~ /[[:space:]]/ {
      value = $NF
      gsub(/[Dd]/, "E", value)
      count++
      magnitude = value + 0.0
      if (magnitude < 0.0) magnitude = -magnitude
      if (value !~ /^[-+]?[0-9.]+[EeDd][-+]?[0-9]+$/ || scale + 0 <= 0 || magnitude >= tolerance * scale) {
        printf "TEST FAILED: %s = %s (tolerance %s x scale %s)\n", label, value, tolerance, scale
        exit 1
      }
    }
    END {
      if (!count) {
        printf "TEST FAILED: missing %s\n", label
        exit 1
      }
    }' "$PROT"
}

check_tensor() {
  label=$1
  pattern="^$label[[:space:]]*[-+0-9]"
  count=$(grep -c "$pattern" "$PROT" || true)
  test "$count" = "3" || {
    echo "TEST FAILED: expected three rows for $label, found $count" >&2
    exit 1
  }
  grep "$pattern" "$PROT" | awk -v label="$label" '
    {
      for (column = NF - 2; column <= NF; column++) {
        value = $column
        numeric = value + 0.0
        if (value !~ /^[-+]?[0-9.]+[EeDd][-+]?[0-9]+$/ || numeric != numeric) {
          printf "TEST FAILED: non-finite %s component %s\n", label, value
          exit 1
        }
      }
    }'
}

test "$(value_for 'NUMBER OF K-POINTS')" = "8" || {
  echo "TEST FAILED: expected eight k-points" >&2
  exit 1
}

check_abs "SKALA POSITIVE TAU CHECK" 1.0e-10
check_abs "PS GRID MINUS TRACE" 1.0e-10
# The CPU iterative orthogonalizer stops at max|S-I| < 1e-8; its occupied
# trace is bounded by that tolerance times the occupation sum.
check_abs "TRACE MINUS OCCUPATIONS" 1.0e-8 "$(value_for 'OCCUPATION ELECTRONS')"
check_abs "GRAD ADJOINT DIFFERENCE" 1.0e-10
check_abs "TAU OPERATOR DIFFERENCE" 1.0e-10
check_abs "ONE-CENTER MATRIX DIFFERENCE" 1.0e-10
check_abs "TAU ANGULAR-RADIAL ERROR" 1.0e-10
# A cold one-step integration test is not an SCF convergence test. Require
# finite residuals and the PAW metric identities, not a small cold residual.
check_abs "SKALA OCCUPIED RESIDUAL RMS" 1.0e99
check_abs "SKALA OCCUPIED RESIDUAL MAX" 1.0e99
check_abs "SKALA OCCUPATION COMMUTATOR MAX" 1.0e99
check_abs "SKALA SCF OVERLAP ERROR" 1.0e-8
check_abs "SKALA HAMILTONIAN HERMITICITY" 1.0e-10
check_tensor "SKALA MODEL D E / D STRAIN"
check_tensor "SKALA CORE D E / D STRAIN"
check_tensor "SKALA TAU D E / D STRAIN"
check_tensor "TOTAL D E / D STRAIN"

value_for "FULL ELECTRONIC OPERATOR CONTRACTION" >/dev/null
grep -q "ANALYTIC STRESS[[:space:]]*AVAILABLE" "$PROT"
grep -q "PROGRAM FINISHED" "$PROT"
echo "TEST PASSED"
