#!/usr/bin/env python3
import argparse
import collections
import csv
import glob
import os
import sys


def profile_files(path):
    matches = glob.glob(path)
    if not matches:
        matches = [path]
    files = []
    for match in matches:
        if os.path.isdir(match):
            files.extend(
                glob.glob(os.path.join(match, "**", "*_profile*.csv"), recursive=True)
            )
        elif os.path.isfile(match):
            files.append(match)
    return sorted(set(files))


def suite_name(root, path):
    if root and os.path.isdir(root):
        return os.path.basename(os.path.normpath(root))
    parts = os.path.normpath(path).split(os.sep)
    if len(parts) >= 4 and parts[-2].startswith("rep"):
        return parts[-4]
    return ""


def case_and_repeat(path):
    parts = os.path.normpath(path).split(os.sep)
    if len(parts) >= 3 and parts[-2].startswith("rep"):
        return parts[-3], parts[-2]
    base = os.path.basename(path)
    case = base.split("_profile", 1)[0]
    return case, ""


def collect(paths):
    rows = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "gbyte": 0.0, "files": set()}
    )
    for arg in paths:
        root = arg if os.path.isdir(arg) else ""
        for path in profile_files(arg):
            suite = suite_name(root, path)
            case, repeat = case_and_repeat(path)
            with open(path, newline="") as handle:
                for row in csv.DictReader(handle):
                    op = row.get("op", "")
                    if not op.startswith("ACC_COPY"):
                        continue
                    gbyte = float(row.get("gbyte") or 0.0)
                    if gbyte <= 0.0:
                        continue
                    key = (suite, case, repeat, op)
                    data = rows[key]
                    data["calls"] += int(row.get("calls") or 0)
                    data["seconds"] += float(row.get("total_seconds") or 0.0)
                    data["gbyte"] += gbyte
                    data["files"].add(path)
    return rows


def aggregate(rows, include_repeat):
    merged = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "gbyte": 0.0, "files": set()}
    )
    for (suite, case, repeat, op), data in rows.items():
        key = (suite, case, repeat if include_repeat else "", op)
        item = merged[key]
        item["calls"] += data["calls"]
        item["seconds"] += data["seconds"]
        item["gbyte"] += data["gbyte"]
        item["files"].update(data["files"])
    return merged


def selected_rows(rows, top, per_case, min_gb):
    items = [
        (key, data)
        for key, data in rows.items()
        if data["gbyte"] >= min_gb
    ]
    if not per_case:
        return sorted(items, key=lambda item: -item[1]["gbyte"])[:top]

    grouped = collections.defaultdict(list)
    for key, data in items:
        suite, case, repeat, _op = key
        grouped[(suite, case, repeat)].append((key, data))
    selected = []
    for group_key in sorted(grouped):
        group = sorted(grouped[group_key], key=lambda item: -item[1]["gbyte"])
        selected.extend(group[:top])
    return selected


def print_tsv(rows):
    print("suite\tcase\trepeat\top\tcalls\tcopy_gb\tseconds\tfiles")
    for (suite, case, repeat, op), data in rows:
        print(
            "{}\t{}\t{}\t{}\t{}\t{:.9g}\t{:.9g}\t{}".format(
                suite,
                case,
                repeat,
                op,
                data["calls"],
                data["gbyte"],
                data["seconds"],
                len(data["files"]),
            )
        )


def print_markdown(rows):
    print("| suite | case | repeat | op | calls | copy_gb | seconds | files |")
    print("| --- | --- | --- | --- | ---: | ---: | ---: | ---: |")
    for (suite, case, repeat, op), data in rows:
        print(
            "| {} | {} | {} | {} | {} | {:.4f} | {:.4f} | {} |".format(
                suite,
                case,
                repeat,
                op,
                data["calls"],
                data["gbyte"],
                data["seconds"],
                len(data["files"]),
            )
        )


def main(argv):
    parser = argparse.ArgumentParser(
        description="Summarize ACC_COPY profile rows from CP-PAW profile CSV files."
    )
    parser.add_argument("paths", nargs="+", help="profile CSV, glob, or run root")
    parser.add_argument("--top", type=int, default=10, help="rows per table/group")
    parser.add_argument(
        "--per-case",
        action="store_true",
        help="emit the top rows separately for each suite/case/repeat group",
    )
    parser.add_argument(
        "--include-repeat",
        action="store_true",
        help="keep repeats separate instead of aggregating them",
    )
    parser.add_argument("--min-gb", type=float, default=0.0)
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args(argv[1:])

    rows = aggregate(collect(args.paths), args.include_repeat)
    chosen = selected_rows(rows, args.top, args.per_case, args.min_gb)
    if not chosen:
        print("No ACC_COPY rows found.", file=sys.stderr)
        return 1
    if args.markdown:
        print_markdown(chosen)
    else:
        print_tsv(chosen)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
