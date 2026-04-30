#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "${HERE}/../../.." && pwd)

cd "${ROOT}"

bash -n paw_install
bash -n src/Buildtools/defaultparmfile
bash -n src/Buildtools/paw_build.sh
bash -n src/Buildtools/paw_fcflags.sh
bash -n src/Buildtools/paw_srclist.sh
bash -n tests/profile/si64/run_benchmark.sh
bash -n tests/profile/si64/run_followup.sh
bash -n tests/profile/si64/run_nsys.sh
bash -n tests/profile/si64/run_overnight.sh

python3 -m py_compile \
  tests/profile/si64/profile_summary.py \
  tests/profile/si64/benchmark_summary.py

test -f tests/profile/si64/si64.cntl
test -f tests/profile/si64/si64.strc
test -f tests/profile/si64/si64_bands.cntl

git diff --check
