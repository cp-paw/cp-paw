#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
if [ "$#" -gt 0 ]; then
  "${MAKE:-make}" -C "$1" -f Makefile -f "$root/tests/unittests/scan/driver.mk" scan-test scan-backend-test
else
  work=$(mktemp -d "${TMPDIR:-/tmp}/scan-kernel.XXXXXX")
  trap 'rm -rf "$work"' EXIT HUP INT TERM
  "${FC:-gfortran}" -O2 -fdollar-ok -fimplicit-none -Wall -Wextra \
    -J"$work" -o "$work/test_scan.x" "$root/src/paw_scan.f90" \
    "$root/tests/unittests/scan/test_scan.f90" "$root/tests/unittests/scan/reference_scan.f90" \
    "$root/tests/unittests/scan/block_scan.f90"
  "$work/test_scan.x"
fi
