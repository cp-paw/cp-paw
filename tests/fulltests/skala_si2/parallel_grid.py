#!/usr/bin/env python3
"""Compare legacy and cooperative Skala grids from one fixed Si2 restart."""

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
from types import SimpleNamespace
import subprocess
import time

from cache_parity import digest
import force_mode_parity as modes
import force_stationary_fd as fd


def counters(text, label, steps, ranks):
    rows = [list(map(int, line[len(label):].split())) for line in text.splitlines()
            if line.startswith(label + " ")]
    if len(rows) != steps or any(len(row) != ranks or min(row) < 0 for row in rows):
        raise ValueError(f"Invalid rank accounting: {label}")
    return rows


def compare(reference, current, mode):
    differences = {}
    if reference["scalars"].keys() != current["scalars"].keys():
        raise ValueError("Changed scalar diagnostics")
    for label in reference["scalars"]:
        if len(reference["scalars"][label]) != len(current["scalars"][label]):
            raise ValueError(f"Changed diagnostic count: {label}")
        differences[label] = max(abs(a-b) for a, b in
                                 zip(reference["scalars"][label], current["scalars"][label]))
        limit = 2.e-5 if "L2" in label or "OPERATOR" in label else 1.e-7
        if not math.isfinite(differences[label]) or differences[label] > limit:
            raise ValueError(f"Changed {label}: {differences[label]}")
    a, b = reference["energies"], current["energies"]
    if len(a) != len(b) or [row["step"] for row in a] != [row["step"] for row in b]:
        raise ValueError("Changed evaluated energy steps")
    differences["energy"] = max(abs(x["energy"]-y["energy"]) for x, y in zip(a, b))
    differences["waves_and_lambda"] = modes.wave_difference(reference["waves"], current["waves"])
    if differences["energy"] > 1.e-7 or differences["waves_and_lambda"] > 1.e-7:
        raise ValueError(f"Changed energy or electronic restart: {differences}")
    if mode != "electronic":
        if len(reference["forces"]) != len(current["forces"]) or any(
                len(a["forces"]) != len(b["forces"]) for a, b in
                zip(reference["forces"], current["forces"])):
            raise ValueError("Changed force layout")
        differences["forces"] = max(abs(x-y) for u, v in zip(reference["forces"], current["forces"])
                                    for a, b in zip(u["forces"], v["forces"]) for x, y in zip(a, b))
        if differences["forces"] > 2.e-6:
            raise ValueError(f"Changed forces: {differences['forces']}")
    if mode == "stress":
        a, b = reference["stress"], current["stress"]
        if len(a) != len(b):
            raise ValueError("Changed stress layout")
        differences["stress"] = max(abs(x-y) for u, v in zip(a, b) for x, y in zip(u, v))
        if differences["stress"] > 2.e-6:
            raise ValueError(f"Changed stress: {differences['stress']}")
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--accelerated-grid", action="store_true")
    parser.add_argument("--root-weight", type=float, default=1.)
    parser.add_argument("--unpacked", action="store_true")
    parser.add_argument("--auto-balance", action="store_true")
    parser.add_argument("--sample-batch", action="store_true")
    parser.add_argument("--sample-gpu", action="store_true")
    parser.add_argument("--tensor-cache-mib", type=int, default=64)
    parser.add_argument("--legacy-uncached", action="store_true")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument("--modes", nargs="+", choices=("electronic", "forces", "stress"),
                        default=["electronic", "forces", "stress"])
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--steps", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--radial", type=int, default=64)
    parser.add_argument("--angular", type=int, default=11)
    parser.add_argument("--cache-mib", type=int, default=256)
    args = parser.parse_args()
    if not math.isfinite(args.root_weight) or not .25 <= args.root_weight <= 16.:
        parser.error("Root weight must lie in 0.25..16")
    if args.auto_balance and (args.device != "CUDA" or args.steps < 4):
        parser.error("Auto balance needs central CUDA inference and at least four steps")
    if args.tensor_cache_mib < 0 or (args.sample_gpu and not args.sample_batch):
        parser.error("Need a nonnegative tensor-cache budget and batched sampling for sample GPU")
    if min(*args.ranks, args.steps, args.repeats, args.radial, args.angular) < 1 or args.cache_mib < 0:
        parser.error("Need positive ranks, step and grid counts, and nonnegative cache budget")
    if args.accelerated_grid and (args.device != "CUDA" or args.gpu_mode == "off" or args.cache_mib == 0):
        parser.error("Accelerated grid checks need CUDA, enabled PAW GPU libraries and a source cache")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    hashes = {key: digest(path) for key, path in inputs.items()}
    env = {key: value for key, value in os.environ.items() if not key.startswith("CPPAW_")}
    env.update(CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_DETERMINISTIC="1", CPPAW_SKALA_SCF_DETAIL="1",
               CPPAW_ACCEL_PROFILE="1", CPPAW_ACCEL_PROFILE_FILE="parallel",
               CPPAW_SKALA_SOURCE_CACHE_MB=str(args.cache_mib), CPPAW_SKALA_PARTITION_CACHE_MB="256",
               CPPAW_SKALA_TORCH_THREADS="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    env.update(CPPAW_SKALA_POINT_ROOT_WEIGHT=str(args.root_weight),
               CPPAW_SKALA_PACKED_EXCHANGE="0" if args.unpacked else "1",
               CPPAW_SKALA_TENSOR_CACHE_MB=str(args.tensor_cache_mib),
               CPPAW_SKALA_SAMPLE_BATCH="1" if args.sample_batch else "0",
               CPPAW_SKALA_SAMPLE_ACC="1" if args.sample_gpu else "0")
    if args.auto_balance:
        env["CPPAW_SKALA_POINT_BALANCE"] = "1"
    if args.accelerated_grid:
        env.update(CPPAW_SKALA_SOURCE_FWD_ACC="1", CPPAW_SKALA_SOURCE_FWD_ACC_MIN_ROWS="1",
                   CPPAW_SKALA_SOURCE_BACK_ACC="1", CPPAW_SKALA_SOURCE_BACK_ACC_MIN_ROWS="1",
                   CPPAW_SKALA_SOURCE_FWD_RESIDENT_MB="256",
                   CPPAW_SKALA_GRID_BACK_ACC="1", CPPAW_SKALA_GRID_BACK_ACC_MIN_POINTS="1")
    fd.write_json(root / "provenance.json", dict(inputs={key: str(path) for key, path in inputs.items()},
        sha256=hashes, arguments={key: str(value) if isinstance(value, Path) else value
                                 for key, value in vars(args).items()},
        environment={key: value for key, value in env.items() if key.startswith("CPPAW_")},
        scope="implementation parity and short timing probes, not physical convergence"))
    geometry = fd.geometry(inputs["restart"].read_bytes())
    results = []
    for mode in args.modes:
        settings = SimpleNamespace(dt=1.e-6, block_steps=args.steps, device=args.device,
            radial_points=args.radial, lebedev_exactness=args.angular, cutoff=20., density_dual=2.,
            mass=25., mass_g2=0.3166286988823056, friction=0.4,
            orthogonality_tolerance=1.e-12, stress=mode == "stress")
        template = fd.control(settings, electronic_only=mode == "electronic")
        for ranks in args.ranks:
            reference = None
            for repeat in range(args.repeats):
                for parallel in (False, True):
                    name = f"{mode}-r{ranks}-{'points' if parallel else 'legacy'}-{repeat+1}"
                    work = root / name
                    work.mkdir()
                    control = template.replace("!SKALA", f"!SKALA PARALLELGRID={'T' if parallel else 'F'}")
                    (work / "si2.cntl").write_text(control)
                    for key, filename in (("structure", "si2.strc"), ("restart", "si2.rstrt")):
                        shutil.copy2(inputs[key], work / filename)
                    (work / "model.fun").symlink_to(inputs["model"])
                    command = [args.mpiexec, "-np", str(ranks), str(inputs["executable"]), "si2.cntl"]
                    print("Running " + name, flush=True)
                    start = time.monotonic()
                    current_env = dict(env)
                    if not parallel and args.legacy_uncached:
                        current_env.update(CPPAW_SKALA_TENSOR_CACHE_MB="0", CPPAW_SKALA_SAMPLE_BATCH="0",
                                           CPPAW_SKALA_SAMPLE_ACC="0", CPPAW_SKALA_PACKED_EXCHANGE="0")
                    with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
                        subprocess.run(command, cwd=work, env=current_env, stdout=out, stderr=err,
                                       timeout=600, check=True)
                    elapsed = time.monotonic()-start
                    text = (work / "si2.prot").read_text()
                    if "PROGRAM FINISHED" not in text:
                        raise ValueError(f"{name}: missing normal termination")
                    restart = (work / "si2.rstrt").read_bytes()
                    fd.check_geometry(geometry, fd.geometry(restart))
                    rows = counters(text, "GRID POINTS PER MPI RANK", args.steps, ranks)
                    calls = counters(text, "MODEL CALLS PER MPI RANK", args.steps, ranks)
                    gpu_rows = {}
                    if args.sample_batch and (parallel or not args.legacy_uncached):
                        sampled = counters(text, "SMOOTH SAMPLE BATCH ROWS", args.steps, 1)
                        if [row[0] for row in sampled] != [sum(row) for row in rows]:
                            raise ValueError("Incomplete batched smooth-grid sampling")
                        if args.sample_gpu:
                            sampled = counters(text, "SMOOTH SAMPLE GPU ROWS", args.steps, 1)
                            if [row[0] for row in sampled] != [row[0] for row in rows]:
                                raise ValueError("GPU smooth-grid sampling missed root-owned points")
                    if args.accelerated_grid:
                        for label in ("SOURCE GPU FORWARD ROWS", "SOURCE GPU REVERSE ROWS"):
                            gpu_rows[label] = [row[0] for row in counters(text, label, args.steps, 1)]
                            if min(gpu_rows[label]) == 0:
                                raise ValueError(f"No actual GPU coverage: {label}")
                        reused = counters(text, "SOURCE GPU REUSED ROWS", args.steps, 1)
                        warm = reused[3:] if args.auto_balance else reused[1:]
                        if args.steps > 1 and min(row[0] for row in warm) == 0:
                            raise ValueError("No warm GPU geometry residency")
                        gpu_rows["SOURCE GPU REUSED ROWS"] = [row[0] for row in reused]
                        copies = 0
                        for path in work.glob("parallel*.csv"):
                            with path.open() as stream:
                                copies += sum(int(row["calls"]) for row in csv.DictReader(stream)
                                    if row["op"] == "ACC_COPY_SKALA_GRID_BACK_OUT")
                        if copies < args.steps:
                            raise ValueError("No device-resident grid-backprojection output")
                    for row, count in zip(rows, calls):
                        if sum(count) != geometry["natom"] or sum(row) < 1:
                            raise ValueError("Missing or duplicated atom-block inference")
                        if args.device == "CUDA" and count[1:] != [0]*(ranks-1):
                            raise ValueError("CUDA inference escaped the central rank")
                        if parallel and sum(row) >= ranks*geometry["natom"] and min(row) == 0:
                            raise ValueError("Cooperative grid left a nonempty point range idle")
                    forces = fd.force_records(text, geometry["natom"]) if mode != "electronic" else []
                    current = dict(scalars=modes.scalars(text, args.steps),
                        energies=forces if forces else modes.energy_records(text),
                        forces=forces, waves=modes.wave_records(restart),
                        stress=[tuple(fd.numbers(line[len("TOTAL D E / D STRAIN"):], 3))
                                for line in text.splitlines()
                                if line.startswith("TOTAL D E / D STRAIN ")])
                    if mode == "stress" and len(current["stress"]) != 3*args.steps:
                        raise ValueError("Missing stress tensor")
                    differences = compare(reference, current, mode) if reference is not None else {}
                    if reference is None:
                        reference = current
                    results.append(dict(name=name, mode=mode, ranks=ranks, parallel=parallel,
                        elapsed_s=elapsed, rows=rows, model_calls=calls, gpu_rows=gpu_rows,
                        differences=differences))
                    fd.write_json(root / "results.json", dict(passed=False, completed=results))
    if any(digest(path) != hashes[key] for key, path in inputs.items()):
        raise ValueError("An executable or scientific input changed during the test")
    fd.write_json(root / "results.json", dict(passed=True, completed=results))
    print("PASS: complete-block inference, rank ownership, energies, restarts and requested derivatives")


if __name__ == "__main__":
    main()
