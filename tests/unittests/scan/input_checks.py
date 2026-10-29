#!/usr/bin/env python3
"""Reject incomplete native SCAN operators before constructing PAW setups."""

import argparse
from contextlib import nullcontext
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("The MPI rank count must be positive")
    command = [str(args.executable.resolve(strict=True)), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    fixture = Path(__file__).resolve().parents[2] / "fulltests/si2"
    template = (fixture / "si2.cntl").read_text()
    cases = [
        ("two-functionals", "!SCAN !END !SKALA MODEL='absent.fun' !END", "", "MUTUALLY EXCLUSIVE"),
        ("model", "!SCAN MODEL='absent.fun' !END", "", "DOES NOT USE A MODEL"),
        ("device", "!SCAN DEVICE='CUDA' !END", "", "DOES NOT USE A MODEL"),
        ("device-index", "!SCAN DEVICEINDEX=0 !END", "", "DOES NOT USE A TORCH DEVICE INDEX"),
        ("apply", "!SCAN APPLY=F !END", "", "REQUIRES APPLY=T"),
        ("smooth", "!SCAN APPLYSMOOTH=F !END", "", "COMPLETE SMOOTH OPERATOR"),
        ("tau", "!SCAN APPLYTAU=F !END", "", "REQUIRES THE TAU OPERATOR"),
        ("one-center", "!SCAN APPLYONECENTER=F !END", "", "COMPLETE PAW OPERATOR"),
        ("distributed-model", "!SCAN DISTRIBUTEGPU=T !END", "", "NO DISTRIBUTED TORCH MODEL"),
        ("dispersion", "!SCAN !END", "VDW=T", "SCAN D3 PARAMETRIZATION IS NOT IMPLEMENTED"),
    ]
    env = dict(os.environ, CPPAW_GPU_MODE="off", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if args.output is not None:
        args.output.mkdir(parents=True, exist_ok=False)
    context = (nullcontext(args.output) if args.output is not None else
               tempfile.TemporaryDirectory(prefix="paw-scan-input-"))
    with context as temporary:
        root = Path(temporary)
        for name, scan, extra, message in cases:
            work = root / name
            work.mkdir()
            for filename in ("si2.strc", "stp.cntl"):
                shutil.copy2(fixture / filename, work / filename)
            text, count = re.subn(r"!DFT\s+TYPE=10\s+!END", f"!DFT TYPE=10 {extra} {scan} !END", template)
            if count != 1:
                raise ValueError("The conventional Si2 fixture changed")
            (work / "si2.cntl").write_text(text)
            run = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=60)
            (work / "stdout.log").write_text(run.stdout)
            (work / "stderr.log").write_text(run.stderr)
            output = run.stdout + run.stderr
            if (work / "si2.prot").exists():
                output += (work / "si2.prot").read_text()
            if message not in output or "PROGRAM FINISHED" in output:
                raise RuntimeError(f"{name}: missing native SCAN rejection {message!r}\n{output[-4000:]}")
            print(f"{name}: rejected as required", flush=True)


if __name__ == "__main__":
    main()
