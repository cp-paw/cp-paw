#!/usr/bin/env python3
"""Compare fixed-cell stress with the heavy-cell reference at matched states."""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

from cache_parity import digest
import force_mode_parity as modes
import force_stationary_fd as fd
import stationarity


def stress_records(text):
    expected = [row["step"] for row in stationarity.records(text)]
    result, current, step = [], None, None
    for line in text.splitlines():
        if line.startswith(stationarity.HEADER):
            if current is not None or (step is not None and
                                       (not result or result[-1]["step"] != step)):
                raise ValueError("Missing or incomplete stress report")
            step = int(line[len(stationarity.HEADER):].strip())
        elif line == "SKALA TOTAL STRESS DIAGNOSTIC":
            if current is not None or step is None or (result and result[-1]["step"] == step):
                raise ValueError("Duplicate or unscoped stress report")
            current = {"step": step, "tensor": []}
        elif line.startswith("TOTAL D E / D STRAIN"):
            if current is None:
                raise ValueError("Unscoped stress tensor row")
            current["tensor"].append(fd.numbers(line[len("TOTAL D E / D STRAIN"):], 3))
            if len(current["tensor"]) == 3:
                result.append(current)
                current = None
    if current is not None or [row["step"] for row in result] != expected:
        raise ValueError("Stress reports do not match electronic steps")
    return result


def stress_difference(first, second):
    if not first or len(first) != len(second):
        raise ValueError("Changed or empty stress trace")
    maximum = 0.
    for a, b in zip(first, second):
        if a["step"] != b["step"] or len(a["tensor"]) != 3 or len(b["tensor"]) != 3:
            raise ValueError("Changed stress step or tensor shape")
        for x, y in zip(a["tensor"], b["tensor"]):
            if len(x) != 3 or len(y) != 3 or not all(math.isfinite(v) for v in (*x, *y)):
                raise ValueError("Malformed or nonfinite stress tensor")
            maximum = max(maximum, *(abs(u-v) for u, v in zip(x, y)))
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--source-cache-mib", type=int, default=0)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--steps", type=int, default=3)
    args = parser.parse_args()
    if args.source_cache_mib < 0 or args.mpi_ranks < 1:
        parser.error("Need a nonnegative cache budget and positive rank count")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("Tolerance must be positive and finite")
    if not math.isfinite(args.dt) or args.dt <= 0 or args.steps < 2:
        parser.error("Need a positive finite timestep and at least two steps")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_DETERMINISTIC="1",
               CPPAW_SKALA_SCF_DETAIL="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               CPPAW_SKALA_SOURCE_CACHE_MB=str(args.source_cache_mib),
               CPPAW_SKALA_PARTITION_CACHE_MB="256")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    settings = SimpleNamespace(dt=args.dt, block_steps=args.steps, device=args.device,
                               radial_points=96, lebedev_exactness=17, cutoff=20.,
                               mass=25., mass_g2=0.3166286988823056, friction=0.4,
                               orthogonality_tolerance=1e-12)
    template = fd.control(settings)
    provenance = {key: {"path": str(path), "sha256": digest(path)} for key, path in inputs.items()}
    provenance["driver_sha256"] = digest(Path(__file__))
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                                for key, value in vars(args).items()}
    provenance["environment"] = {k: v for k, v in env.items() if k.startswith("CPPAW_") or
                                 k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")}
    fd.write_json(args.output / "provenance.json", provenance)
    command = [str(inputs["executable"]), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    before = fd.geometry(inputs["restart"].read_bytes())
    result, waves, stresses = {}, {}, {}
    cases = (("no-stress", "FORCE=T STRESS=F", False),
             ("fixed-cell", "FORCE=T STRESS=T", False),
             ("required-forces", "FORCE=F STRESS=T", False),
             ("heavy-cell", "FORCE=F STRESS=F", True))
    for name, flags, moving in cases:
        work = args.output / name
        work.mkdir()
        control = template.replace("!PSIDYN", f"!PSIDYN {flags}")
        if moving:
            control = control.replace("!CELL MOVE=F", "!CELL MOVE=T")
        (work / "si2.cntl").write_text(control)
        (work / "model.fun").symlink_to(inputs["model"])
        shutil.copy2(inputs["structure"], work / "si2.strc")
        shutil.copy2(inputs["restart"], work / "si2.rstrt")
        fd.write_json(work / "inputs.json", {key: digest(work / key)
                      for key in ("si2.cntl", "si2.strc", "si2.rstrt", "model.fun")})
        print(f"Running {name}: {work}", flush=True)
        with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
            subprocess.run(command, cwd=work, env=env, stdout=out, stderr=err,
                           timeout=1800, check=True)
        text = (work / "si2.prot").read_text()
        trace = stationarity.records(text)
        stationarity.validate(trace)
        if len(trace) != settings.block_steps:
            raise ValueError("Unexpected number of electronic steps")
        result[name] = dict(trace=trace, bands=stationarity.band_records(text),
                            energies=fd.force_records(text, before["natom"]),
                            scalars=modes.scalars(text, len(trace)))
        if name == "no-stress":
            if "SKALA TOTAL STRESS DIAGNOSTIC" in text or "TOTAL D E / D STRAIN" in text:
                raise ValueError("Unrequested stress report")
        else:
            stresses[name] = stress_records(text)
        restart = (work / "si2.rstrt").read_bytes()
        fd.check_geometry(before, fd.geometry(restart))
        waves[name] = modes.wave_records(restart)
    differences = {}
    for name, _, _ in cases[1:]:
        row = modes.compare(result["no-stress"], result[name], args.tolerance)
        row["wave_and_lambda"] = modes.wave_difference(waves["no-stress"], waves[name],
                                                      cell_tolerance=1e-12 if name == "heavy-cell" else 0.)
        row["total_force"] = max(abs(x-y)
            for a, b in zip(result["no-stress"]["energies"], result[name]["energies"])
            for u, v in zip(a["forces"], b["forces"]) for x, y in zip(u, v))
        row["strain_derivative"] = stress_difference(stresses["heavy-cell"], stresses[name])
        if any(value > args.tolerance for value in row.values()):
            raise ValueError(f"{name}: stress-mode mismatch {row}")
        differences[name] = row
    for key, path in inputs.items():
        if digest(path) != provenance[key]["sha256"]:
            raise ValueError(f"Changed input during verification: {key}")
    summary = dict(passed=True, differences=differences, tolerance=args.tolerance,
                   scope="fixed-cell stress implementation parity, not stationary strain accuracy",
                   stress_trace=stresses)
    fd.write_json(args.output / "results.json", summary)
    print(json.dumps({k:v for k,v in summary.items() if k != "stress_trace"}, indent=2))


if __name__ == "__main__":
    main()
