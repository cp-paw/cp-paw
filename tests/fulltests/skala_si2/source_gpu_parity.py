#!/usr/bin/env python3
"""Compare bounded GPU source reverse with CPU evaluation at matched states."""

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
from stress_mode_parity import stress_difference, stress_records


def gpu_rows(text, steps):
    label = "SOURCE GPU REVERSE ROWS"
    rows = [fd.numbers(line[len(label):], 1)[0] for line in text.splitlines()
            if line.startswith(label + " ")]
    if len(rows) != steps or any(x < 0 or not x.is_integer() for x in rows):
        raise ValueError("Missing or invalid GPU row accounting")
    return [int(x) for x in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--modes", nargs="+", choices=("electronic", "force", "stress"),
                        default=["electronic", "force", "stress"])
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    if args.mpi_ranks < 1 or len(set(args.modes)) != len(args.modes):
        parser.error("Need positive rank count and distinct calculation modes")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("Tolerance must be positive and finite")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_DETERMINISTIC="1",
               CPPAW_SKALA_SCF_DETAIL="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               CPPAW_SKALA_SOURCE_BACK_ACC_MIN_ROWS="1", CPPAW_SKALA_PARTITION_CACHE_MB="256")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    settings = SimpleNamespace(dt=0.01, block_steps=2, device=args.device,
                               radial_points=96, lebedev_exactness=17, cutoff=20.,
                               mass=25., mass_g2=0.3166286988823056, friction=0.4,
                               orthogonality_tolerance=1e-12)
    template = fd.control(settings)
    provenance = {key: {"path": str(path), "sha256": digest(path)} for key, path in inputs.items()}
    provenance["driver_sha256"] = digest(Path(__file__))
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                                for key, value in vars(args).items()}
    provenance["environment"] = {key: value for key, value in env.items()
        if key.startswith("CPPAW_") or key in
        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")}
    fd.write_json(args.output / "provenance.json", provenance)
    command = [str(inputs["executable"]), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    before = fd.geometry(inputs["restart"].read_bytes())
    variants = (("off", "0", "256", "4096"), ("gpu", "1", "256", "4096"),
                ("bounded", "1", "1", "4096"), ("zero-budget", "1", "0", "4096"),
                ("no-cache", "1", "256", "0"))
    summaries = {}
    for mode in args.modes:
        flags = dict(electronic="FORCE=F STRESS=F", force="FORCE=T STRESS=F",
                     stress="FORCE=T STRESS=T")[mode]
        result, waves, stresses, coverage = {}, {}, {}, {}
        for name, enabled, budget, cache in variants:
            work = args.output / mode / name
            work.mkdir(parents=True)
            (work / "si2.cntl").write_text(template.replace("!PSIDYN", f"!PSIDYN {flags}"))
            (work / "model.fun").symlink_to(inputs["model"])
            shutil.copy2(inputs["structure"], work / "si2.strc")
            shutil.copy2(inputs["restart"], work / "si2.rstrt")
            fd.write_json(work / "inputs.json", {key: digest(work / key)
                          for key in ("si2.cntl", "si2.strc", "si2.rstrt", "model.fun")})
            current_env = dict(env, CPPAW_SKALA_SOURCE_BACK_ACC=enabled,
                               CPPAW_SKALA_SOURCE_BACK_ACC_MB=budget,
                               CPPAW_SKALA_SOURCE_CACHE_MB=cache)
            fd.write_json(work / "environment.json", {key: value for key, value in current_env.items()
                          if key.startswith("CPPAW_")})
            print(f"Running {mode}/{name}: {work}", flush=True)
            with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
                subprocess.run(command, cwd=work, env=current_env, stdout=out, stderr=err,
                               timeout=1800, check=True)
            text = (work / "si2.prot").read_text()
            trace = stationarity.records(text)
            stationarity.validate(trace)
            if len(trace) != settings.block_steps:
                raise ValueError("Unexpected number of electronic steps")
            result[name] = dict(trace=trace, bands=stationarity.band_records(text),
                energies=modes.energy_records(text) if mode == "electronic" else
                         fd.force_records(text, before["natom"]),
                scalars=modes.scalars(text, len(trace)))
            coverage[name] = gpu_rows(text, len(trace))
            if mode == "stress":
                stresses[name] = stress_records(text)
            elif "SKALA TOTAL STRESS DIAGNOSTIC" in text:
                raise ValueError("Unrequested stress calculation")
            restart = (work / "si2.rstrt").read_bytes()
            fd.check_geometry(before, fd.geometry(restart))
            waves[name] = modes.wave_records(restart)
        if any(x != 0 for name in ("off", "zero-budget", "no-cache") for x in coverage[name]):
            raise ValueError("Disabled GPU path processed rows")
        if any(not 0 < small < full for small, full in zip(coverage["bounded"], coverage["gpu"])):
            raise ValueError("Expected both full and partial device coverage")
        differences = {}
        for name, _, _, _ in variants[1:]:
            row = modes.compare(result["off"], result[name], args.tolerance)
            row["wave_and_lambda"] = modes.wave_difference(waves["off"], waves[name])
            if mode != "electronic":
                row["total_force"] = max(abs(x-y)
                    for a, b in zip(result["off"]["energies"], result[name]["energies"])
                    for u, v in zip(a["forces"], b["forces"]) for x, y in zip(u, v))
            if mode == "stress":
                row["strain_derivative"] = stress_difference(stresses["off"], stresses[name])
            if any(value > args.tolerance for value in row.values()):
                raise ValueError(f"{mode}/{name}: source GPU mismatch {row}")
            differences[name] = row
        summaries[mode] = dict(differences=differences, gpu_rows=coverage)
        print(json.dumps({mode: summaries[mode]}), flush=True)
    for key, path in inputs.items():
        if digest(path) != provenance[key]["sha256"]:
            raise ValueError(f"Changed input during verification: {key}")
    fd.write_json(args.output / "results.json", dict(passed=True, modes=summaries,
        tolerance=args.tolerance, scope="source reverse implementation parity, not physical convergence"))


if __name__ == "__main__":
    main()
