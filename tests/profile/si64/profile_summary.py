#!/usr/bin/env python3
import collections
import csv
import glob
import sys


def category(op):
    if op.startswith("MPI_ALLTOALL"):
        return "MPI alltoall"
    if op.startswith("FFT") or op.startswith("PW_FFT"):
        return "FFT"
    if "GEMM" in op or "HERK" in op or "SYRK" in op:
        return "BLAS-like"
    return "other"


def main(argv):
    files = []
    for arg in argv[1:]:
        files.extend(glob.glob(arg))
    files = sorted(set(files))
    if not files:
        print("No profile CSV files found.", file=sys.stderr)
        return 1

    per_op = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "gflop": 0.0, "gbyte": 0.0}
    )
    per_shape = collections.defaultdict(
        lambda: {"calls": 0, "seconds": 0.0, "max": 0.0, "gflop": 0.0}
    )

    for path in files:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                op = row["op"]
                calls = int(row["calls"])
                seconds = float(row["total_seconds"])
                gflop = float(row["gflop"])
                gbyte = float(row["gbyte"])
                key = (op, row["n1"], row["n2"], row["n3"], row["n4"])

                op_data = per_op[op]
                op_data["calls"] += calls
                op_data["seconds"] += seconds
                op_data["gflop"] += gflop
                op_data["gbyte"] += gbyte

                shape_data = per_shape[key]
                shape_data["calls"] += calls
                shape_data["seconds"] += seconds
                shape_data["max"] = max(shape_data["max"], float(row["max_seconds"]))
                shape_data["gflop"] += gflop

    by_category = collections.defaultdict(float)
    for op, data in per_op.items():
        by_category[category(op)] += data["seconds"]
    total = sum(by_category.values())

    print("Profile files: {}".format(len(files)))
    print("Instrumented rank-seconds: {:.6f}".format(total))
    print("")
    print("Category summary")
    for name, seconds in sorted(by_category.items(), key=lambda item: -item[1]):
        percent = 100.0 * seconds / total if total else 0.0
        print("  {:<14s} {:10.4f}s {:6.2f}%".format(name, seconds, percent))

    print("")
    print("Top operations")
    for op, data in sorted(per_op.items(), key=lambda item: -item[1]["seconds"])[:15]:
        seconds = data["seconds"]
        rate = data["gflop"] / seconds if seconds else 0.0
        bw = data["gbyte"] / seconds if seconds else 0.0
        print(
            "  {:<24s} calls={:8d} seconds={:9.4f} GF/s={:8.2f} GB/s={:8.2f}".format(
                op, data["calls"], seconds, rate, bw
            )
        )

    print("")
    print("Top shapes")
    for key, data in sorted(per_shape.items(), key=lambda item: -item[1]["seconds"])[:15]:
        op, n1, n2, n3, n4 = key
        seconds = data["seconds"]
        rate = data["gflop"] / seconds if seconds else 0.0
        print(
            "  {:<24s} n=({},{},{},{}) calls={:7d} seconds={:9.4f} GF/s={:8.2f}".format(
                op, n1, n2, n3, n4, data["calls"], seconds, rate
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
