#!/usr/bin/env python3
"""Verify optional force suppression and required dynamics using conventional Si2."""

import argparse
from contextlib import nullcontext
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "fulltests/skala_si2"))
from force_mode_parity import wave_records, wave_difference, trajectory_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--output", type=Path, help="Retain all inputs and outputs for inspection")
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("MPI rank count must be positive")
    command = [str(args.executable.resolve(strict=True)), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    fixture = Path(__file__).resolve().parents[2] / "fulltests/si2"
    template = (fixture / "si2.cntl").read_text().replace("NSTEP=180", "NSTEP=2")
    template = template.replace("!CONTROL\n", "!CONTROL\n !ANALYSE !TRA FORCE=T E=T !END !END\n")
    env = dict(os.environ, CPPAW_GPU_MODE="off", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    if args.output is not None:
        args.output.mkdir(parents=True, exist_ok=False)
    context = (nullcontext(args.output) if args.output is not None else
               tempfile.TemporaryDirectory(prefix="paw-force-selection-"))
    with context as temporary:
        root = Path(temporary)
        cases = [("seed", "", ""), ("default", "", ""),
                 ("explicit", "FORCE=T", ""), ("omitted", "FORCE=F", ""),
                 ("atoms-on", "FORCE=T", "!RDYN FRIC=0.0 !END"),
                 ("atoms-off", "FORCE=F", "!RDYN FRIC=0.0 !END"),
                 ("cell-on", "FORCE=T", "!CELL MOVE=T FRIC=0.0 M=1.E30 !END"),
                 ("cell-off", "FORCE=F STRESS=F", "!CELL MOVE=T FRIC=0.0 M=1.E30 !END"),
                 ("stress-on", "FORCE=T STRESS=T", "!CELL MOVE=F FRIC=0.0 M=1.E30 !END"),
                 ("stress-off", "FORCE=F STRESS=T", "!CELL MOVE=F FRIC=0.0 M=1.E30 !END")]
        waves = {}
        for name, settings, dynamics in cases:
            work = root / name
            work.mkdir()
            for filename in ("si2.strc", "stp.cntl"):
                shutil.copy2(fixture / filename, work / filename)
            text = template.replace("!PSIDYN", f"{dynamics}\n !PSIDYN {settings}")
            if name != "seed":
                text = text.replace("START=t", "START=f")
                shutil.copy2(root / "seed/si2.rstrt", work / "si2.rstrt")
            (work / "si2.cntl").write_text(text)
            run = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=300)
            (work / "stdout.log").write_text(run.stdout)
            (work / "stderr.log").write_text(run.stderr)
            output = run.stdout + run.stderr
            if (work / "si2.prot").exists():
                output += (work / "si2.prot").read_text()
            if run.returncode != 0 or "PROGRAM FINISHED" not in output:
                raise RuntimeError(f"{name} failed:\n{output[-4000:]}")
            trajectory = work / "si2_f.tra"
            has_forces = trajectory.exists() and trajectory.stat().st_size > 0
            if has_forces != (name != "omitted"):
                raise RuntimeError(f"{name}: incorrect force-trajectory availability")
            waves[name] = wave_records((work / "si2.rstrt").read_bytes())
            print(f"{name}: passed", flush=True)
        for a, b in (("default", "explicit"), ("explicit", "omitted"),
                     ("atoms-on", "atoms-off"), ("cell-on", "cell-off"),
                     ("stress-on", "stress-off"), ("explicit", "stress-on")):
            difference = wave_difference(waves[a], waves[b],
                                         cell_tolerance=1e-12 if a.startswith("cell") else 0.)
            if difference > 1e-10:
                raise RuntimeError(f"{a}/{b}: electronic restart difference {difference}")
            energy = trajectory_difference((root / a / "si2_e.tra").read_bytes(),
                                           (root / b / "si2_e.tra").read_bytes())
            if energy > 1e-10:
                raise RuntimeError(f"{a}/{b}: energy trajectory difference {energy}")
            force = 0.
            if a.startswith(("atoms", "cell", "stress")):
                force = trajectory_difference((root / a / "si2_f.tra").read_bytes(),
                                              (root / b / "si2_f.tra").read_bytes())
                if force > 1e-10:
                    raise RuntimeError(f"{a}/{b}: required force difference {force}")
            print(f"{a}/{b}: wave/Lambda {difference:.3e}, energy {energy:.3e}, force {force:.3e}", flush=True)


if __name__ == "__main__":
    main()
