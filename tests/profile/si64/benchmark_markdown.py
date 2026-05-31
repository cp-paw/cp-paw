#!/usr/bin/env python3
import csv
import os
import re
import sys


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
    print("| suite | case | ranks | ok | wall_s | rank_s | gap_s | coverage_% | paw_s | blas_s | lapack_s | fft_s | mpi_s | pw_trace_s | phase_s | phase_gap_s | copy_gb | energy | energy_delta |")
    print("| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in rows:
        print(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                suite_label(row),
                row.get("case", ""),
                row.get("ranks", ""),
                row.get("ok", ""),
                number(row.get("wall_s")),
                number(row.get("rank_s")),
                number(row.get("gap_s")),
                number(row.get("coverage_pct")),
                number(row.get("paw_s")),
                number(row.get("blas_s")),
                number(row.get("lapack_s")),
                number(row.get("fft_s")),
                number(row.get("mpi_s")),
                number(row.get("pw_trace_s")),
                number(row.get("phase_s")),
                number(row.get("phase_gap_s")),
                number(row.get("copy_gb")),
                number(row.get("energy"), 6),
                number(row.get("energy_delta"), 6),
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
