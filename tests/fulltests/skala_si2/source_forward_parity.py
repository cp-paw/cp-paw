#!/usr/bin/env python3
"""Check source-forward offload alone and with source reverse at matched states."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import force_mode_parity as modes
import force_stationary_fd as fd
import source_gpu_parity as reverse
import stationarity
from stress_mode_parity import stress_difference, stress_records


def coverage_rows(text, steps, label):
    rows = [fd.numbers(line[len(label):], 1)[0] for line in text.splitlines()
            if line.startswith(label + " ")]
    if len(rows) != steps or any(not math.isfinite(x) or x < 0 or not x.is_integer() for x in rows):
        raise ValueError("Missing or invalid " + label)
    return [int(x) for x in rows]


def forward_rows(text, steps):
    return coverage_rows(text, steps, "SOURCE GPU FORWARD ROWS")


def reused_rows(text, steps):
    return coverage_rows(text, steps, "SOURCE GPU REUSED ROWS")


def check_reuse(forward, reused, expected):
    if len(forward) != len(reused) or not reused or reused[0] != 0:
        raise ValueError("Invalid first-use residency accounting")
    if any(not 0 <= a <= b for a, b in zip(reused, forward)):
        raise ValueError("Resident coverage exceeds forward coverage")
    if expected == "full":
        if len(reused) < 2 or any(a != b or a == 0 for a, b in zip(reused[1:], forward[1:])):
            raise ValueError("Expected complete warm residency")
    elif expected == "bounded":
        if len(reused) < 2 or any(not 0 < a < b for a, b in zip(reused[1:], forward[1:])):
            raise ValueError("Expected positive, pool-limited warm residency")
    elif expected != "off" or any(reused):
        raise ValueError("Unexpected residency")


def check_resident_profile(work, ranks, forward, reused):
    paths = sorted(work.glob("source_forward_profile*.csv"))
    if len(paths) != ranks:
        raise ValueError("Unexpected profile rank count")
    names = ("ACC_COPY_SKALA_FWD_GEOM", "ACC_COPY_SKALA_FWD_INPUT",
             "SKALA_SOURCE_FWD_DEVICE", "ACC_COPY_SKALA_FWD_OUT",
             "ACC_PRESENT_SKALA_FWD_GEOM")
    events = {name: dict(calls=0, rows=0, bytes=0.) for name in names}
    for path in paths:
        with path.open() as stream:
            for row in csv.DictReader(stream):
                if row["op"] not in events:
                    continue
                calls, np, n, ne, nc = (int(row[key]) for key in ("calls", "n1", "n2", "n3", "n4"))
                size = float(row["gbyte"]) * 1e9
                if min(calls, np, n, ne) <= 0 or nc != 4 or not math.isfinite(size) or size < 0:
                    raise ValueError("Invalid source residency profile record")
                if row["op"] == names[0]:
                    expected = calls * (4 * (np + 1) + 8 * ne * (8 * n + 5))
                    if not math.isclose(size, expected, rel_tol=5e-9, abs_tol=0.5):
                        raise ValueError("Incorrect geometry upload payload")
                if row["op"] == names[4] and size != 0:
                    raise ValueError("Resident geometry reported duplicate transfer bytes")
                if row["op"] == names[3] and not math.isclose(size, 80 * np * calls,
                        rel_tol=5e-9, abs_tol=0.5):
                    raise ValueError("Incorrect forward field output payload")
                event = events[row["op"]]
                event["calls"] += calls
                event["rows"] += np * calls
                event["bytes"] += size
    total, hits = sum(forward), sum(reused)
    if events[names[0]]["rows"] != total - hits:
        raise ValueError("Geometry uploads disagree with residency counters")
    if any(events[name]["rows"] != total for name in names[1:4]):
        raise ValueError("Electronic inputs, kernel or output missing from a call")
    if not hits <= events[names[4]]["rows"] <= total:
        raise ValueError("Invalid present-geometry accounting")
    return events


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--residency", action="store_true",
                        help="Check repeated forward calls with full, bounded and disabled GPU pools")
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("Need a positive rank count")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    env = {key: value for key, value in os.environ.items() if not key.startswith("CPPAW_")}
    env.update(CPPAW_GPU_MODE="off", CPPAW_SKALA_DETERMINISTIC="1",
               CPPAW_SKALA_SCF_DETAIL="1", CPPAW_SKALA_GRID_BACK_ACC="1",
               CPPAW_SKALA_GRID_BACK_ACC_MIN_POINTS="1", CPPAW_ACCEL_PROFILE="1",
               CPPAW_ACCEL_PROFILE_FILE="source_forward_profile",
               CPPAW_SKALA_SOURCE_FWD_ACC_MIN_ROWS="1", CPPAW_SKALA_SOURCE_BACK_ACC_MIN_ROWS="1",
               CPPAW_SKALA_SOURCE_BACK_ACC_MB="4096", CPPAW_SKALA_PARTITION_CACHE_MB="256",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    settings = SimpleNamespace(dt=5., block_steps=2, device=args.device,
        radial_points=96, lebedev_exactness=17, cutoff=20., mass=25.,
        mass_g2=0.3166286988823056, friction=0.1, orthogonality_tolerance=1e-12)
    meta = dict(inputs={key: dict(path=str(path), sha256=fd.digest(path)) for key, path in inputs.items()},
        driver_sha256=fd.digest(Path(__file__)), tolerance=1e-10,
        settings=vars(settings), mpi_ranks=args.mpi_ranks, residency=args.residency,
        environment={key: value for key, value in env.items() if key.startswith("CPPAW_") or
                     key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")},
        scope="Matched-state forward/reverse backend parity, not physical convergence or timings")
    fd.write_json(output / "provenance.json", meta)
    command = [str(inputs["executable"]), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    before = fd.geometry(inputs["restart"].read_bytes())
    variants = (("off", "0", "4096", "4096", "0"),
                ("forward", "1", "4096", "4096", "0"),
                ("bounded", "1", "1", "4096", "0"),
                ("zero-budget", "1", "0", "4096", "0"),
                ("no-cache", "1", "4096", "0", "0"),
                ("reverse", "0", "4096", "4096", "1"),
                ("both", "1", "4096", "4096", "1"))
    variants = tuple((*item, "0", "off") for item in variants)
    if args.residency:
        variants = (("off", "0", "4096", "4096", "0", "4096", "off"),
                    ("forward", "1", "4096", "4096", "0", "0", "off"),
                    ("resident", "1", "4096", "4096", "0", "4096", "full"),
                    ("resident-pool-bounded", "1", "4096", "4096", "0", "256", "bounded"),
                    ("bounded", "1", "1", "4096", "0", "4096", "full"),
                    ("zero-budget", "1", "0", "4096", "0", "4096", "off"),
                    ("no-cache", "1", "4096", "0", "0", "4096", "off"),
                    ("both", "1", "4096", "4096", "1", "4096", "full"))
    results = {}
    for mode, flags in (("electronic", "FORCE=F STRESS=F"), ("force", "FORCE=T STRESS=F"),
                        ("stress", "FORCE=T STRESS=T")):
        baseline = None
        coverage = {}
        for name, enabled, budget, cache, back, pool, expect_reuse in variants:
            work = output / mode / name
            work.mkdir(parents=True)
            (work / "si2.cntl").write_text(fd.control(settings).replace("!PSIDYN", f"!PSIDYN {flags}"))
            (work / "model.fun").symlink_to(inputs["model"])
            shutil.copy2(inputs["structure"], work / "si2.strc")
            shutil.copy2(inputs["restart"], work / "si2.rstrt")
            fd.write_json(work / "inputs.json", {key: fd.digest(work / key)
                for key in ("si2.cntl", "si2.strc", "si2.rstrt", "model.fun")})
            selected = dict(env, CPPAW_SKALA_SOURCE_FWD_ACC=enabled,
                CPPAW_SKALA_SOURCE_FWD_ACC_MB=budget, CPPAW_SKALA_SOURCE_CACHE_MB=cache,
                CPPAW_SKALA_SOURCE_BACK_ACC=back, CPPAW_SKALA_SOURCE_FWD_RESIDENT_MB=pool)
            fd.write_json(work / "environment.json", {key: value for key, value in selected.items()
                                                       if key.startswith("CPPAW_")})
            print(f"Running source-forward parity {mode}/{name}", flush=True)
            with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
                subprocess.run(command, cwd=work, env=selected, stdout=out, stderr=err,
                               check=True, timeout=1800)
            text = (work / "si2.prot").read_text()
            trace = stationarity.records(text)
            stationarity.validate(trace)
            if len(trace) != settings.block_steps:
                raise ValueError("Unexpected number of evaluations")
            data = (work / "si2.rstrt").read_bytes()
            fd.check_geometry(before, fd.geometry(data))
            actual = dict(trace=trace, bands=stationarity.band_records(text),
                energies=modes.energy_records(text) if mode == "electronic" else
                         fd.force_records(text, before["natom"]),
                scalars=modes.scalars(text, len(trace)), waves=modes.wave_records(data),
                stress=stress_records(text) if mode == "stress" else [])
            if mode != "stress" and "SKALA TOTAL STRESS DIAGNOSTIC" in text:
                raise ValueError("Unrequested stress evaluation")
            fwd = forward_rows(text, len(trace))
            rev = reverse.gpu_rows(text, len(trace))
            if any(x != 0 for x in fwd) and (enabled == "0" or budget == "0" or cache == "0"):
                raise ValueError("Disabled forward path processed rows")
            if any(x != 0 for x in rev) != (back == "1"):
                raise ValueError("Unexpected reverse coverage")
            coverage[name] = fwd
            if baseline is None:
                baseline = actual
            differences = modes.compare(baseline, actual, meta["tolerance"])
            differences["wave_and_lambda"] = modes.wave_difference(baseline["waves"], actual["waves"])
            if mode != "electronic":
                differences["total_force"] = max(abs(x-y)
                    for a, b in zip(baseline["energies"], actual["energies"])
                    for u, v in zip(a["forces"], b["forces"]) for x, y in zip(u, v))
            if mode == "stress":
                differences["strain_derivative"] = stress_difference(baseline["stress"], actual["stress"])
            row = dict(differences=differences, forward_rows=fwd, reverse_rows=rev,
                       passed=all(math.isfinite(x) and x <= meta["tolerance"] for x in differences.values()))
            if args.residency:
                row["reused_rows"] = reused_rows(text, len(trace))
                check_reuse(fwd, row["reused_rows"], expect_reuse)
                row["profile"] = check_resident_profile(work, args.mpi_ranks, fwd, row["reused_rows"])
            results[f"{mode}/{name}"] = row
            fd.write_json(output / "results.json", results)
            print(json.dumps(row), flush=True)
            if not row["passed"]:
                raise ValueError("Unchanged source-forward parity gate failed")
        if any(not 0 < a < b for a, b in zip(coverage["bounded"], coverage["forward"])):
            raise ValueError("Expected positive and budget-limited forward coverage")
        if coverage["forward"] != coverage["both"]:
            raise ValueError("Reverse changed forward coverage")
        if args.residency and any(coverage[name] != coverage["forward"]
                for name in ("resident", "resident-pool-bounded")):
            raise ValueError("Residency changed forward coverage")
    for name, path in inputs.items():
        if fd.digest(path) != meta["inputs"][name]["sha256"]:
            raise ValueError("Changed input during verification")
    fd.write_json(output / "complete.json", dict(passed=True, cases=len(results), tolerance=meta["tolerance"]))


if __name__ == "__main__":
    main()
