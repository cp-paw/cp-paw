#!/usr/bin/env python3
import collections
import os
import statistics
import sys

from benchmark_compare import parse_float, read_rows, speedup, suite_label


BASELINES = (
    ("vs_1cpu_gnu", "cpu1", "cpu"),
    ("vs_1cpu_nvhpc", "cpu1", "nvhpc_cpu"),
    ("vs_8cpu_nvhpc", "cpu8", "nvhpc_cpu"),
)


def number(value, digits=2):
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def valid_row(row):
    return row.get("ok") == "yes" and row.get("_wall") is not None


def median_wall(rows):
    values = [row["_wall"] for row in rows if valid_row(row)]
    return statistics.median(values) if values else None


def baseline_walls(rows):
    walls = {}
    groups = {row["_group"] for row in rows}
    for group in groups:
        group_rows = [row for row in rows if row["_group"] == group]
        for label, kind, case in BASELINES:
            candidates = [
                row
                for row in group_rows
                if row["_kind"] == kind and row.get("case") == case
            ]
            walls[(group, label)] = median_wall(candidates)
    return walls


def aggregate(rows):
    groups = collections.defaultdict(list)
    for row in rows:
        key = (
            row["_group"],
            row.get("suite", ""),
            row.get("case", ""),
            row.get("ranks", ""),
        )
        groups[key].append(row)

    baselines = baseline_walls(rows)
    result = []
    for (group, suite, case, ranks), case_rows in groups.items():
        valid = [row for row in case_rows if valid_row(row)]
        walls = [row["_wall"] for row in valid]
        energy_deltas = [
            abs(value)
            for value in (parse_float(row.get("energy_delta")) for row in valid)
            if value is not None
        ]
        median = statistics.median(walls) if walls else None
        result.append(
            {
                "group": group,
                "suite": suite,
                "case": case,
                "ranks": ranks,
                "valid": len(valid),
                "total": len(case_rows),
                "median": median,
                "minimum": min(walls) if walls else None,
                "maximum": max(walls) if walls else None,
                "energy_delta_max": max(energy_deltas) if energy_deltas else None,
                **{
                    label: speedup(baselines.get((group, label)), median)
                    for label, _, _ in BASELINES
                },
            }
        )
    return sorted(
        result,
        key=lambda row: (
            row["group"],
            {"gpu": 0, "cpu1": 1, "cpu8": 2}.get(
                next(
                    (
                        item["_kind"]
                        for item in rows
                        if item["_group"] == row["group"]
                        and item.get("suite", "") == row["suite"]
                        and item.get("case", "") == row["case"]
                        and item.get("ranks", "") == row["ranks"]
                    ),
                    "",
                ),
                3,
            ),
            row["case"],
            row["suite"],
        ),
    )


def main(argv):
    if len(argv) < 2:
        print("usage: benchmark_medians.py benchmark.tsv [...]", file=sys.stderr)
        return 2

    rows = []
    for path in argv[1:]:
        rows.extend(read_rows(path))
    if not rows:
        return 1

    print("### CP-PAW Benchmark Medians")
    print()
    sources = ", ".join(f"`{os.path.basename(path)}`" for path in argv[1:])
    print(f"Sources: {sources}")
    print()
    print(
        "| suite | case | ranks | valid | median_s | min_s | max_s | "
        "vs_1cpu_gnu | vs_1cpu_nvhpc | vs_8cpu_nvhpc | energy_delta_max |"
    )
    print(
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    )
    for row in aggregate(rows):
        print(
            "| {} | {} | {} | {}/{} | {} | {} | {} | {} | {} | {} | {} |".format(
                suite_label({"suite": row["suite"]}),
                row["case"],
                row["ranks"],
                row["valid"],
                row["total"],
                number(row["median"]),
                number(row["minimum"]),
                number(row["maximum"]),
                number(row["vs_1cpu_gnu"]),
                number(row["vs_1cpu_nvhpc"]),
                number(row["vs_8cpu_nvhpc"]),
                number(row["energy_delta_max"], 8),
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
