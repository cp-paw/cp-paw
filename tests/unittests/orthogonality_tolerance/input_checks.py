#!/usr/bin/env python3
"""Exercise PSIDYN tolerance parsing in isolated conventional Si2 calculations."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("MPI rank count must be positive")
    executable = args.executable.resolve(strict=True)
    command = [str(executable), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    fixture = Path(__file__).resolve().parents[2] / "fulltests/si2"
    control = (fixture / "si2.cntl").read_text().replace("NSTEP=180", "NSTEP=1")
    cases = [
        ("default", "", None),
        ("explicit-default", "SAFEORTHO=T ORTHOTOL=1.E-8", None),
        ("tight", "SAFEORTHO=T ORTHOTOL=1.E-12", None),
        ("alternative-default", "SAFEORTHO=F", None),
        ("incompatible", "SAFEORTHO=F ORTHOTOL=1.E-12", "ORTHOTOL REQUIRES SAFEORTHO=T"),
        ("out-of-range", "SAFEORTHO=T ORTHOTOL=1.E-7", "ORTHOTOL MUST LIE BETWEEN"),
    ]
    env = dict(os.environ, CPPAW_GPU_MODE="off", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="paw-orthogonality-") as temporary:
        for name, settings, error in cases:
            work = Path(temporary) / name
            work.mkdir()
            for filename in ("si2.strc", "stp.cntl"):
                shutil.copy2(fixture / filename, work / filename)
            (work / "si2.cntl").write_text(control.replace("!PSIDYN", f"!PSIDYN {settings}"))
            run = subprocess.run(command, cwd=work, env=env,
                                 capture_output=True, text=True, timeout=180)
            protocol = work / "si2.prot"
            text = run.stdout + run.stderr + (protocol.read_text() if protocol.exists() else "")
            if error is None:
                if run.returncode != 0 or "PROGRAM FINISHED" not in text:
                    raise RuntimeError(f"{name} failed:\n{text[-4000:]}")
            elif run.returncode == 0 or error not in text:
                raise RuntimeError(f"{name} did not reject invalid input:\n{text[-4000:]}")
            print(f"{name}: passed", flush=True)


if __name__ == "__main__":
    main()
