#!/bin/bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
build=${1:?Usage: run.sh path/to/bin/Build_profile}
make_tool=${MAKE:-make}
"$make_tool" -C "$build" -f Makefile -f "$root/tests/unittests/skala_reconstruction/driver.mk" skala-reconstruction-test
