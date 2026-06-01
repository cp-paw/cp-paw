#!/usr/bin/env python3
import csv
import glob
import os
import re
import sys


def copy_bucket(op):
    if "OFFDEN" in op:
        return "copy_offden_gb"
    if "DENMAT" in op:
        return "copy_denmat_gb"
    if any(token in op for token in ("PROPSI", "PRO_CACHE", "THIS_PROJ")):
        return "copy_proj_gb"
    if any(
        token in op
        for token in (
            "_PSI0",
            "_PSI1",
            "_PSI2",
            "_PSI_IN",
            "_PSI_OUT",
            "_PSIM",
            "_OPSI",
            "_HPSI",
        )
    ):
        return "copy_wave_gb"
    if any(token in op for token in ("PROJ", "ADDPRO")):
        return "copy_proj_gb"
    return None


def profile_totals(run_dir):
    totals = {
        "instrumented": 0.0,
        "blas": 0.0,
        "lapack": 0.0,
        "fft": 0.0,
        "mpi": 0.0,
        "paw": 0.0,
        "pw_trace": 0.0,
        "pw_gtor": 0.0,
        "pw_rtog": 0.0,
        "phase": 0.0,
        "setup": 0.0,
        "copy_gb": 0.0,
        "copy_wave_gb": 0.0,
        "copy_proj_gb": 0.0,
        "copy_offden_gb": 0.0,
        "copy_denmat_gb": 0.0,
    }
    for path in glob.glob(os.path.join(run_dir, "*_profile*.csv")):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                op = row["op"]
                seconds = float(row["total_seconds"])
                gbyte = float(row["gbyte"])
                if op.startswith("ACC_COPY"):
                    totals["copy_gb"] += gbyte
                    bucket = copy_bucket(op)
                    if bucket:
                        totals[bucket] += gbyte
                    continue
                if op.startswith("ACC_SETUP"):
                    totals["setup"] += seconds
                    continue
                if op.startswith("PHASE_"):
                    totals["phase"] += seconds
                    continue
                if op.startswith("PW_") and not op.startswith("PW_FFT"):
                    totals["pw_trace"] += seconds
                    if op.startswith("PW_GTOR_"):
                        totals["pw_gtor"] += seconds
                    elif op.startswith("PW_RTOG_"):
                        totals["pw_rtog"] += seconds
                    continue
                if op.startswith("PAW_"):
                    totals["paw"] += seconds
                    continue
                totals["instrumented"] += seconds
                if op.startswith("MPI_ALLTOALL"):
                    totals["mpi"] += seconds
                elif op.startswith("FFT") or op.startswith("PW_FFT") or op.startswith("CUFFT"):
                    totals["fft"] += seconds
                elif op.startswith("LAPACK") or op.startswith("CUSOLVER"):
                    totals["lapack"] += seconds
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


def energy_check(energy, env):
    expected = env.get("expected_energy") or env.get("EXPECTED_ENERGY")
    if not expected:
        return None, None
    if energy is None:
        return None, False
    try:
        expected_value = float(expected)
        tolerance = float(env.get("energy_tol") or env.get("ENERGY_TOL") or "1e-5")
    except ValueError:
        return None, False
    delta = abs(energy - expected_value)
    return delta, delta <= tolerance


def wall_rank_time(wall, ranks):
    if wall is None:
        return None
    try:
        return wall * int(ranks)
    except (TypeError, ValueError):
        return None


def main(argv):
    root = argv[1] if len(argv) > 1 else "."
    rows = []
    for run_dir in sorted(glob.glob(os.path.join(root, "*", "rep*"))):
        case = os.path.basename(os.path.dirname(run_dir))
        repeat = os.path.basename(run_dir)
        totals = profile_totals(run_dir)
        env = run_env(run_dir)
        wall = wall_time(run_dir)
        wall_rank = wall_rank_time(wall, env.get("ranks"))
        energy = final_energy(run_dir)
        energy_delta, energy_is_ok = energy_check(energy, env)
        ok = run_ok(run_dir) and (energy_is_ok is not False)
        gap = None
        coverage = None
        phase_gap = None
        if wall_rank is not None:
            gap = wall_rank - totals["instrumented"]
            phase_gap = wall_rank - totals["phase"]
            if wall_rank > 0.0:
                coverage = 100.0 * totals["instrumented"] / wall_rank
        rows.append(
            {
                "case": case,
                "repeat": repeat,
                "nsteps": env.get("nsteps") or env.get("nstps"),
                "ranks": env.get("ranks"),
                "ok": "yes" if ok else "no",
                "wall_s": wall,
                "wall_rank_s": wall_rank,
                "rank_s": totals["instrumented"],
                "gap_s": gap,
                "coverage_pct": coverage,
                "paw_s": totals["paw"],
                "blas_s": totals["blas"],
                "lapack_s": totals["lapack"],
                "fft_s": totals["fft"],
                "mpi_s": totals["mpi"],
                "pw_trace_s": totals["pw_trace"],
                "pw_gtor_s": totals["pw_gtor"],
                "pw_rtog_s": totals["pw_rtog"],
                "phase_s": totals["phase"],
                "phase_gap_s": phase_gap,
                "setup_s": totals["setup"],
                "copy_gb": totals["copy_gb"],
                "copy_wave_gb": totals["copy_wave_gb"],
                "copy_proj_gb": totals["copy_proj_gb"],
                "copy_offden_gb": totals["copy_offden_gb"],
                "copy_denmat_gb": totals["copy_denmat_gb"],
                "energy": energy,
                "energy_delta": energy_delta,
                "energy_ok": (
                    "" if energy_is_ok is None else ("yes" if energy_is_ok else "no")
                ),
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
        "wall_rank_s",
        "rank_s",
        "gap_s",
        "coverage_pct",
        "paw_s",
        "blas_s",
        "lapack_s",
        "fft_s",
        "mpi_s",
        "pw_trace_s",
        "pw_gtor_s",
        "pw_rtog_s",
        "phase_s",
        "phase_gap_s",
        "setup_s",
        "copy_gb",
        "copy_wave_gb",
        "copy_proj_gb",
        "copy_offden_gb",
        "copy_denmat_gb",
        "energy",
        "energy_delta",
        "energy_ok",
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
