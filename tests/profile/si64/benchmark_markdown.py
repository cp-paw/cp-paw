#!/usr/bin/env python3
import csv
import os
import re
import sys


def parse_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def number(value, digits=2):
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except ValueError:
        return str(value)


def suite_label(row):
    suite = row.get("suite", "")
    if not suite:
        return ""
    match = re.search(r"empty(\d+).*nstep(\d+)", suite)
    if match:
        return f"empty={match.group(1)}, nsteps={match.group(2)}"
    return suite


def read_rows(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if "suite" not in row:
            row["suite"] = ""
        if row.get("transfer_gb") in (None, ""):
            copy_gb = parse_float(row.get("copy_gb"))
            update_gb = parse_float(row.get("update_gb"))
            if copy_gb is not None or update_gb is not None:
                row["transfer_gb"] = str((copy_gb or 0.0) + (update_gb or 0.0))
    return rows


def main(argv):
    if len(argv) != 2:
        print("usage: benchmark_markdown.py benchmark.tsv", file=sys.stderr)
        return 2
    path = argv[1]
    rows = read_rows(path)
    if not rows:
        return 1

    print(f"### CP-PAW Benchmark Summary")
    print()
    print(f"Source: `{os.path.basename(path)}`")
    print()
    columns = [
        ("suite", None),
        ("case", None),
        ("ranks", None),
        ("ok", None),
        ("wall_s", 2),
        ("rank_s", 2),
        ("gap_s", 2),
        ("coverage_pct", 2),
        ("paw_s", 2),
        ("skala_s", 2),
        ("skala_partition_s", 2),
        ("skala_atom_grid_s", 2),
        ("skala_onecenter_s", 2),
        ("skala_model_s", 2),
        ("skala_onecenter_adjoint_s", 2),
        ("skala_grid_back_s", 2),
        ("blas_s", 2),
        ("lapack_s", 2),
        ("fft_s", 2),
        ("fft_kernel_s", 2),
        ("mpi_s", 2),
        ("pw_trace_s", 2),
        ("pw_gtor_s", 2),
        ("pw_rtog_s", 2),
        ("phase_s", 2),
        ("phase_gap_s", 2),
        ("transfer_gb", 2),
        ("copy_gb", 2),
        ("copy_wave_gb", 2),
        ("copy_proj_gb", 2),
        ("copy_offden_gb", 2),
        ("copy_denmat_gb", 2),
        ("update_gb", 2),
        ("update_wave_gb", 2),
        ("update_proj_gb", 2),
        ("update_offden_gb", 2),
        ("update_denmat_gb", 2),
        ("energy", 6),
        ("static_total_energy", 6),
        ("model_xc_energy", 6),
        ("hybrid_grid_rows", 0),
        ("partition_classes", 0),
        ("energy_delta", 6),
    ]
    labels = ["coverage_%" if key == "coverage_pct" else key for key, _ in columns]
    print("| " + " | ".join(labels) + " |")
    left_aligned = {"suite", "case", "ok"}
    print(
        "| "
        + " | ".join("---" if key in left_aligned else "---:" for key, _ in columns)
        + " |"
    )
    for row in rows:
        values = []
        for key, digits in columns:
            value = suite_label(row) if key == "suite" else row.get(key, "")
            values.append(str(value) if digits is None else number(value, digits))
        print("| " + " | ".join(values) + " |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
