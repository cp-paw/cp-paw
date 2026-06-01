#!/usr/bin/env python3
import csv
import math
import os
import re
import sys


CPU_BASELINE_CASES = ("nvhpc_cpu", "cpu")


def parse_float(value):
    if value is None or value == "":
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    if math.isnan(result):
        return None
    return result


def number(value, digits=2):
    if value is None:
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def speedup(base_wall, wall):
    if base_wall is None or wall is None or wall <= 0.0:
        return None
    return base_wall / wall


def suite_group_kind(suite):
    if not suite:
        return "", ""
    match = re.match(r"^gpu_acc_(.*)_1rank$", suite)
    if match:
        return f"gpu_acc_{match.group(1)}", ""
    match = re.match(r"^cpu_ref_(.*)_\d+ranks?$", suite)
    if match:
        return f"gpu_acc_{match.group(1)}", "cpu8"
    match = re.match(r"^cusolver_(.*)_1rank$", suite)
    if match:
        return f"cusolver_{match.group(1)}", ""
    match = re.match(r"^cusolver_cpu_ref_(.*)_\d+ranks?$", suite)
    if match:
        return f"cusolver_{match.group(1)}", "cpu8"
    match = re.match(r"^gpu_(\d+)rank$", suite)
    if match:
        return "default", "gpu"
    match = re.match(r"^cpu_1rank$", suite)
    if match:
        return "default", "cpu1"
    match = re.match(r"^cpu_(\d+)rank_ref$", suite)
    if match:
        return "default", "cpu8"
    rank_patterns = (
        (r"^(.*)_(\d+)ranks?_gpu$", "gpu"),
        (r"^(.*)_(\d+)ranks?_cpu_ref$", "cpu8"),
        (r"^(.*)_(\d+)ranks?_cpu$", "cpu"),
    )
    for pattern, kind in rank_patterns:
        match = re.match(pattern, suite)
        if match:
            if kind == "cpu":
                kind = "cpu1" if match.group(2) == "1" else "cpu8"
            return match.group(1), kind
    match = re.match(r"^(.*)-1r$", suite)
    if match:
        return match.group(1), "gpu"
    if suite == "one_rank_gpu":
        return "default", "gpu"
    if suite == "one_rank_cpu":
        return "default", "cpu1"
    if suite == "eight_rank_cpu":
        return "default", "cpu8"
    if "gpu" in suite:
        return suite, "gpu"
    if "cpu_ref" in suite or "eight_rank_cpu" in suite:
        return suite, "cpu8"
    if "cpu" in suite:
        return suite, "cpu1"
    return suite, ""


def case_kind(row):
    group, kind = suite_group_kind(row.get("suite", ""))
    if kind:
        return group, kind
    case = row.get("case", "")
    ranks = row.get("ranks", "")
    if case in CPU_BASELINE_CASES:
        return group, "cpu8" if ranks == "8" else "cpu1"
    if case.startswith(("gpu", "cublas", "cufft", "cusolver")):
        return group, "gpu"
    return group, ""


def read_rows(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row.setdefault("suite", "")
        row["_wall"] = parse_float(row.get("wall_s"))
        row["_group"], row["_kind"] = case_kind(row)
    return rows


def baseline_for(rows, kind):
    candidates = [
        row for row in rows
        if row.get("_kind") == kind
        and row.get("ok", "") == "yes"
        and row.get("case") in CPU_BASELINE_CASES
        and row.get("_wall") is not None
    ]
    if not candidates:
        candidates = [
            row for row in rows
            if row.get("_kind") == kind
            and row.get("ok", "") == "yes"
            and row.get("_wall") is not None
        ]
    if not candidates:
        return None
    return min(candidates, key=lambda row: row["_wall"])


def suite_label(row):
    suite = row.get("suite", "")
    if not suite:
        return ""
    match = re.search(r"empty(\d+).*nstep(\d+)", suite)
    if match:
        return f"empty={match.group(1)}, nsteps={match.group(2)}"
    return suite


def main(argv):
    if len(argv) != 2:
        print("usage: benchmark_compare.py combined_benchmark.tsv", file=sys.stderr)
        return 2
    path = argv[1]
    rows = read_rows(path)
    if not rows:
        return 1

    groups = {}
    for row in rows:
        groups.setdefault(row["_group"], []).append(row)

    print("### CP-PAW Benchmark Comparison")
    print()
    print(f"Source: `{os.path.basename(path)}`")
    print()
    print(
        "| suite | case | ranks | ok | wall_s | vs_1cpu | vs_8cpu | "
        "rank_s | copy_gb | copy_wave_gb | copy_proj_gb | copy_offden_gb | "
        "copy_denmat_gb | energy_delta |"
    )
    print(
        "| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | "
        "---: | ---: | ---: | ---: | ---: |"
    )
    for group in sorted(groups):
        group_rows = groups[group]
        one_cpu = baseline_for(group_rows, "cpu1")
        eight_cpu = baseline_for(group_rows, "cpu8")
        one_cpu_wall = None if one_cpu is None else one_cpu["_wall"]
        eight_cpu_wall = None if eight_cpu is None else eight_cpu["_wall"]
        for row in sorted(
            group_rows,
            key=lambda item: (
                {"gpu": 0, "cpu1": 1, "cpu8": 2}.get(item.get("_kind"), 3),
                item.get("case", ""),
                item.get("ranks", ""),
            ),
        ):
            print(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
                    suite_label(row),
                    row.get("case", ""),
                    row.get("ranks", ""),
                    row.get("ok", ""),
                    number(row.get("_wall")),
                    number(speedup(one_cpu_wall, row.get("_wall"))),
                    number(speedup(eight_cpu_wall, row.get("_wall"))),
                    number(row.get("rank_s")),
                    number(row.get("copy_gb")),
                    number(row.get("copy_wave_gb")),
                    number(row.get("copy_proj_gb")),
                    number(row.get("copy_offden_gb")),
                    number(row.get("copy_denmat_gb")),
                    number(row.get("energy_delta"), 6),
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
