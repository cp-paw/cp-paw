#!/usr/bin/env python3
"""Compare native SCAN energies and complete operators across MPI layouts."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys


def difference(a, b):
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            raise ValueError("Changed result dimensions across MPI layouts")
        return max((difference(x, y) for x, y in zip(a, b)), default=0.)
    if not all(isinstance(x, (int, float)) and math.isfinite(x) for x in (a, b)):
        raise ValueError("Nonfinite or malformed parallel comparison")
    return abs(a-b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("serial-executable", "mpi-executable", "restart", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", default=(2, 3))
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args()
    if (min(args.ranks) < 2 or len(set(args.ranks)) != len(args.ranks)
            or not math.isfinite(args.tolerance) or args.tolerance <= 0):
        parser.error("Need distinct MPI counts greater than one and a positive finite tolerance")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    cases = [(1, args.serial_executable), *((n, args.mpi_executable) for n in args.ranks)]
    results = {}
    for ranks, executable in cases:
        output = args.output / f"ranks-{ranks}"
        command = [sys.executable, str(Path(__file__).with_name("run.py")),
                   "--executable", str(executable.resolve(strict=True)),
                   "--restart", str(args.restart.resolve(strict=True)),
                   "--output", str(output), "--mpi-ranks", str(ranks), "--stress"]
        with (args.output / f"ranks-{ranks}.log").open("w") as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                           timeout=1800, check=True)
        results[ranks] = json.loads((output / "results.json").read_text())
    reference = results[1]
    differences = {}
    for ranks, result in results.items():
        row = {key: difference(reference[key], result[key]) for key in (
            "native_scan_xc_hartree", "reconstructed_xc_hartree", "frozen_core_xc_hartree")}
        for key in ("energy", "forces"):
            row[key] = difference(reference["final"][key], result["final"][key])
        row["strain_derivative"] = difference(reference["stress"][0]["tensor"], result["stress"][0]["tensor"])
        row["operator_norm"] = max(difference(reference["operators"][key], result["operators"][key])
                                   for key in reference["operators"])
        if any(value > args.tolerance for value in row.values()):
            raise ValueError(f"Native SCAN MPI layout {ranks} differs: {row}")
        differences[ranks] = row
    summary = dict(passed=True, tolerance=args.tolerance, differences=differences,
                   scope="matched-state MPI parity, not electronic stationarity or physical accuracy")
    (args.output / "results.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
