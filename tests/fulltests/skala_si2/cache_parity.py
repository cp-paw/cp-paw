#!/usr/bin/env python3
"""Compare two applied-Skala steps with partition caching disabled/enabled."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import stationarity


SCALARS = (*stationarity.LABELS, "MODEL XC ENERGY", "COMPOSITE ELECTRONS",
           "COMPOSITE POSITIVE TAU", "SMOOTH SCALAR OPERATOR L2",
           "SMOOTH TAU OPERATOR L2", "ONE-CENTER MATRIX DIFFERENCE")
LAYOUT = {label: (1, 2) for label in SCALARS}
LAYOUT.update({"TOTAL ENERGY": (1, 2), "SKALA FORCE ATOM": (4, 4),
               "LOCAL PARTITION FORCE ATOM": (4, 4), "TOTAL FORCE ATOM": (4, 4),
               "TOTAL D E / D STRAIN": (3, 6), "SKALA PARTITION STRESS": (3, 6)})
COUNTERS = ("ATOM-GRID ROWS", "PARTITION CACHE HITS", "PARTITION CACHE MISSES",
            "PARTITION CACHE BYTES ALL RANKS")


def diagnostics(text, stress=False):
    data = stationarity.records(text)
    stationarity.validate(data)
    if len(data) != 2:
        raise ValueError("Expected exactly two electronic diagnostic steps")
    layout = {label: shape for label, shape in LAYOUT.items()
              if stress or label not in ("TOTAL D E / D STRAIN", "SKALA PARTITION STRESS")}
    result = {label: [] for label in (*layout, *COUNTERS)}
    force_section = False
    for line in text.splitlines():
        if line == "SKALA TOTAL FORCE DIAGNOSTIC":
            force_section = True
            continue
        if force_section and line.startswith("NET FORCE"):
            force_section = False
        if force_section and re.match(r"^ATOM\s", line):
            line = "TOTAL FORCE ATOM " + line[4:]
        for label in result:
            if not re.match(re.escape(label) + r"\s", line):
                continue
            raw = line[len(label):].split()
            if label == "TOTAL ENERGY" and raw and raw[0] == ":":
                continue
            if not raw or any(not stationarity.NUMBER.fullmatch(x) for x in raw):
                raise ValueError(f"Malformed {label}: {line}")
            values = tuple(float(x.replace("D", "E").replace("d", "e")) for x in raw)
            if not all(math.isfinite(x) for x in values):
                raise ValueError(f"Nonfinite {label}")
            if label in COUNTERS and any(x < 0 or not x.is_integer() for x in values):
                raise ValueError(f"Invalid integer counter {label}")
            width = LAYOUT[label][0] if label in LAYOUT else 1
            if len(values) != width:
                raise ValueError(f"Wrong number of components for {label}")
            result[label].append(values)
    for label, rows in result.items():
        count = LAYOUT[label][1] if label in LAYOUT else 2
        if len(rows) != count:
            raise ValueError(f"Expected {count} rows for {label}, got {len(rows)}")
    return result


def compare(off, cached, tolerance, require_reuse=True, total_force_tolerance=1e-8):
    if any(not math.isfinite(x) or x <= 0 for x in (tolerance, total_force_tolerance)):
        raise ValueError("Tolerance must be positive and finite")
    rows = cached["ATOM-GRID ROWS"]
    if rows != off["ATOM-GRID ROWS"] or rows[0] != rows[1] or rows[0][0] <= 0:
        raise ValueError("Changed or empty grid in cache parity probe")
    zero = [(0.0,), (0.0,)]
    if off["PARTITION CACHE HITS"] != zero or off["PARTITION CACHE MISSES"] != rows:
        raise ValueError("Disabled cache did not recompute every row")
    if off["PARTITION CACHE BYTES ALL RANKS"] != zero:
        raise ValueError("Disabled cache allocated payload")
    hits, misses = cached["PARTITION CACHE HITS"], cached["PARTITION CACHE MISSES"]
    if hits[0] != (0.0,) or any(h[0] + m[0] != n[0] for h, m, n in zip(hits, misses, rows)):
        raise ValueError("Invalid cache hit/miss accounting")
    if require_reuse and (hits[1] != rows[1] or misses[1] != (0.0,)):
        raise ValueError("Second step did not reuse every cached row")
    sizes = cached["PARTITION CACHE BYTES ALL RANKS"]
    if sizes[0] != sizes[1] or sizes[0][0] <= 0:
        raise ValueError("Invalid cache allocation diagnostic")
    differences = {}
    if set(off) != set(cached):
        raise ValueError("Diagnostic layouts differ")
    for label in sorted(off.keys() - set(COUNTERS)):
        delta = max(abs(a - b) for ra, rb in zip(off[label], cached[label])
                    for a, b in zip(ra, rb))
        limit = total_force_tolerance if label == "TOTAL FORCE ATOM" else tolerance
        if delta > limit:
            raise ValueError(f"{label}: cache difference {delta:.6e} > {limit:.6e}")
        differences[label] = delta
    return {"scope": "partition cache parity, not scientific convergence",
            "tolerance": tolerance, "total_force_tolerance": total_force_tolerance,
            "max_abs_difference": differences, "rows_per_step": int(rows[0][0]),
            "cache_hits_per_step": [int(x[0]) for x in hits],
            "cache_misses_per_step": [int(x[0]) for x in misses],
            "cache_bytes_all_ranks": int(sizes[0][0])}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--radial-points", type=int, default=96)
    parser.add_argument("--lebedev-exactness", type=int, default=17)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--total-force-tolerance", type=float, default=1e-8)
    parser.add_argument("--timestep", type=float, default=0.001)
    parser.add_argument("--stress", action="store_true",
                        help="Enable the heavy-cell stress probe; tiny cell changes can invalidate caching")
    args = parser.parse_args()
    if min(args.mpi_ranks, args.radial_points, args.lebedev_exactness) < 1:
        parser.error("Ranks and grid sizes must be positive")
    if not math.isfinite(args.timestep) or args.timestep <= 0:
        parser.error("Timestep must be positive and finite")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    args.output.mkdir(parents=True, exist_ok=False)
    provenance = {key: {"path": str(path), "sha256": digest(path)}
                  for key, path in inputs.items()}
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                               for key, value in vars(args).items()}
    provenance["environment"] = {key: value for key, value in os.environ.items()
                                  if key.startswith("CPPAW_") or key in
                                  ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    template = Path(__file__).with_name("skala_si2.cntl").read_text()
    template = template.replace("NSTEP=1", "NSTEP=2").replace("START=T", "START=F")
    template = template.replace("DT=0.000001", f"DT={args.timestep:.16g}")
    if not args.stress:
        template = template.replace("!CELL MOVE=T", "!CELL MOVE=F")
    template = template.replace("DEVICE='AUTO'", f"DEVICE='{args.device}'")
    template = template.replace("RADIALPOINTS=200", f"RADIALPOINTS={args.radial_points}")
    template = template.replace("LEBEDEVEXACTNESS=53", f"LEBEDEVEXACTNESS={args.lebedev_exactness}")
    template = template.replace("NAME='../si2/stp.cntl'", "NAME='stp.cntl'")
    results = []
    timings = {}
    for name, budget in (("off", "0"), ("cached", "256")):
        work = (args.output / name).resolve()
        work.mkdir()
        (work / "si2.cntl").write_text(template)
        (work / "model.fun").symlink_to(inputs["model"])
        shutil.copy2(inputs["restart"], work / "si2.rstrt")
        shutil.copy2(inputs["structure"], work / "si2.strc")
        shutil.copy2(Path(__file__).parent / "../si2/stp.cntl", work / "stp.cntl")
        env = dict(os.environ, CPPAW_SKALA_PARTITION_CACHE_MB=budget)
        command = [str(inputs["executable"]), "si2.cntl"]
        if args.mpi_ranks > 1:
            command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
        print(f"Running {name}: {work}", flush=True)
        started = time.monotonic()
        with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
            subprocess.run(command, cwd=work, env=env, stdout=out, stderr=err, check=True)
        timings[name] = time.monotonic() - started
        results.append(diagnostics((work / "si2.prot").read_text(), stress=args.stress))
    summary = compare(*results, args.tolerance, require_reuse=not args.stress,
                      total_force_tolerance=args.total_force_tolerance)
    summary["wall_seconds_including_startup_and_first_fill"] = timings
    (args.output / "results.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
