#!/usr/bin/env python3
import csv
import glob
import os
import re
import sys


def profile_totals(run_dir):
    totals = {"instrumented": 0.0, "blas": 0.0, "fft": 0.0, "mpi": 0.0}
    for path in glob.glob(os.path.join(run_dir, "*_profile*.csv")):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                op = row["op"]
                seconds = float(row["total_seconds"])
                totals["instrumented"] += seconds
                if op.startswith("MPI_ALLTOALL"):
                    totals["mpi"] += seconds
                elif op.startswith("FFT") or op.startswith("PW_FFT"):
                    totals["fft"] += seconds
                elif "GEMM" in op or "HERK" in op or "SYRK" in op:
                    totals["blas"] += seconds
    return totals


def wall_time(run_dir):
    path = os.path.join(run_dir, "time.txt")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        for line in handle:
            match = re.match(r"real\s+([0-9.]+)", line)
            if match:
                return float(match.group(1))
    return None


def final_energy(run_dir):
    path = os.path.join(run_dir, "out.log")
    energy = None
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as handle:
        for line in handle:
            if "CONSTANT ENERGY" in line:
                values = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", line)
                if values:
                    energy = float(values[0])
    return energy


def run_env(run_dir):
    values = {}
    path = os.path.join(run_dir, "run.env")
    if not os.path.exists(path):
        return values
    with open(path, errors="replace") as handle:
        for line in handle:
            key, sep, value = line.rstrip("\n").partition("=")
            if sep:
                values[key] = value
    return values


def run_ok(run_dir):
    path = os.path.join(run_dir, "out.log")
    if not os.path.exists(path):
        return False
    with open(path, errors="replace") as handle:
        text = handle.read()
    if "ERROR CODE:" in text:
        return "ERROR CODE:         0" in text
    return "NORMAL STOP" in text


def main(argv):
    root = argv[1] if len(argv) > 1 else "."
    rows = []
    for run_dir in sorted(glob.glob(os.path.join(root, "*", "rep*"))):
        case = os.path.basename(os.path.dirname(run_dir))
        repeat = os.path.basename(run_dir)
        totals = profile_totals(run_dir)
        env = run_env(run_dir)
        rows.append(
            {
                "case": case,
                "repeat": repeat,
                "nsteps": env.get("nsteps") or env.get("nstps"),
                "ranks": env.get("ranks"),
                "ok": "yes" if run_ok(run_dir) else "no",
                "wall_s": wall_time(run_dir),
                "rank_s": totals["instrumented"],
                "blas_s": totals["blas"],
                "fft_s": totals["fft"],
                "mpi_s": totals["mpi"],
                "energy": final_energy(run_dir),
                "env": env.get("env"),
            }
        )

    if not rows:
        print("No benchmark runs found.", file=sys.stderr)
        return 1

    fields = [
        "case",
        "repeat",
        "nsteps",
        "ranks",
        "ok",
        "wall_s",
        "rank_s",
        "blas_s",
        "fft_s",
        "mpi_s",
        "energy",
        "env",
    ]
    print("\t".join(fields))
    for row in rows:
        values = []
        for field in fields:
            value = row[field]
            if isinstance(value, float):
                values.append(f"{value:.9g}")
            elif value is None:
                values.append("")
            else:
                values.append(str(value))
        print("\t".join(values))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
