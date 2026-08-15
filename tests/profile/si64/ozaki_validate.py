#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
import re
import sys


NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def as_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def output_metrics(run_dir):
    metrics = {
        "energy": None,
        "kpoints": None,
        "ortho": [],
        "force_active": None,
        "stress_active": None,
        "forces": {},
        "stress": {},
    }
    out_path = os.path.join(run_dir, "out.log")
    if os.path.exists(out_path):
        with open(out_path, errors="replace") as handle:
            for line in handle:
                if "CONSTANT ENERGY" in line:
                    values = re.findall(NUMBER, line)
                    if values:
                        metrics["energy"] = as_float(values[0])
                match = re.match(r"CPPAW FP64 ORTHO RESIDUAL\s+(\S+)", line)
                if match:
                    metrics["ortho"].append(as_float(match.group(1)))
                match = re.match(r"CPPAW FP64 FORCE ACTIVE\s+([TF])", line)
                if match:
                    metrics["force_active"] = match.group(1) == "T"
                match = re.match(
                    rf"CPPAW FP64 FORCE\s+(\d+)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})",
                    line,
                )
                if match:
                    metrics["forces"][int(match.group(1))] = tuple(
                        as_float(value) for value in match.groups()[1:]
                    )
                match = re.match(r"CPPAW FP64 STRESS ACTIVE\s+([TF])", line)
                if match:
                    metrics["stress_active"] = match.group(1) == "T"
                match = re.match(
                    rf"CPPAW FP64 STRESS\s+(\d+)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})",
                    line,
                )
                if match:
                    metrics["stress"][int(match.group(1))] = tuple(
                        as_float(value) for value in match.groups()[1:]
                    )

    for path in sorted(glob.glob(os.path.join(run_dir, "*.prot"))):
        with open(path, errors="replace") as handle:
            for line in handle:
                if "NUMBER OF K-POINTS" in line:
                    values = re.findall(r"\d+", line)
                    if values:
                        metrics["kpoints"] = int(values[-1])
    return metrics


def ozaki_dispatch(run_dir):
    dispatch = {
        kernel: {"used": 0, "fallback": 0, "unknown": 0, "max_bits": None}
        for kernel in ("DGEMM", "ZGEMM", "ZHERK")
    }
    pattern = re.compile(
        r"^CUBLAS_OZAKI_(DGEMM|ZGEMM|ZHERK)_(USED|FALLBACK|UNKNOWN)$"
    )
    for path in glob.glob(os.path.join(run_dir, "*_profile*.csv")):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                match = pattern.match(row["op"])
                if not match:
                    continue
                kernel, result = match.groups()
                dispatch[kernel][result.lower()] += int(row["calls"])
                if result == "USED":
                    bits = int(row["n1"])
                    current = dispatch[kernel]["max_bits"]
                    dispatch[kernel]["max_bits"] = (
                        bits if current is None else max(current, bits)
                    )
    return dispatch


def max_difference(left, right):
    if left.keys() != right.keys():
        return None
    result = 0.0
    for key in left:
        for a, b in zip(left[key], right[key]):
            result = max(result, abs(a - b))
    return result


def check_limit(label, value, limit, failures):
    if value is None or not math.isfinite(value):
        failures.append(f"{label}: unavailable")
    elif value > limit:
        failures.append(f"{label}: {value:.6e} exceeds {limit:.6e}")


def main(argv):
    parser = argparse.ArgumentParser(
        description="Compare native-FP64 and cuBLAS Ozaki CP-PAW runs."
    )
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--energy-tol", type=float, default=1.0e-8)
    parser.add_argument("--ortho-tol", type=float, default=1.0e-8)
    parser.add_argument("--force-tol", type=float, default=1.0e-7)
    parser.add_argument("--stress-tol", type=float, default=1.0e-7)
    parser.add_argument("--expected-kpoints", type=int)
    parser.add_argument("--require-force", action="store_true")
    parser.add_argument("--require-stress", action="store_true")
    parser.add_argument(
        "--require-ozaki",
        action="append",
        choices=("DGEMM", "ZGEMM", "ZHERK"),
        default=[],
    )
    args = parser.parse_args(argv[1:])

    reference = output_metrics(args.reference)
    candidate = output_metrics(args.candidate)
    dispatch = ozaki_dispatch(args.candidate)
    failures = []

    energy_delta = None
    if reference["energy"] is not None and candidate["energy"] is not None:
        energy_delta = abs(candidate["energy"] - reference["energy"])
    check_limit("energy difference", energy_delta, args.energy_tol, failures)

    if not reference["ortho"] or not candidate["ortho"]:
        failures.append("orthonormality residual: unavailable")
        ortho_max = None
    else:
        ortho_max = max(abs(value) for value in candidate["ortho"])
        check_limit("orthonormality residual", ortho_max, args.ortho_tol, failures)

    if reference["kpoints"] != candidate["kpoints"]:
        failures.append(
            f"k-point count differs: {reference['kpoints']} != {candidate['kpoints']}"
        )
    if (
        args.expected_kpoints is not None
        and candidate["kpoints"] != args.expected_kpoints
    ):
        failures.append(
            f"k-point count {candidate['kpoints']} != {args.expected_kpoints}"
        )

    force_delta = max_difference(reference["forces"], candidate["forces"])
    if args.require_force:
        if not reference["force_active"] or not candidate["force_active"]:
            failures.append("force diagnostic: inactive")
        check_limit("force difference", force_delta, args.force_tol, failures)

    stress_delta = max_difference(reference["stress"], candidate["stress"])
    if args.require_stress:
        if not reference["stress_active"] or not candidate["stress_active"]:
            failures.append("stress diagnostic: inactive")
        check_limit("stress difference", stress_delta, args.stress_tol, failures)

    for kernel in args.require_ozaki:
        if dispatch[kernel]["used"] == 0:
            failures.append(f"{kernel}: no Ozaki dispatch observed")

    print(
        "ozaki_validation "
        f"energy_delta={energy_delta!s} ortho_max={ortho_max!s} "
        f"force_delta={force_delta!s} stress_delta={stress_delta!s} "
        f"kpoints={candidate['kpoints']}"
    )
    for kernel in ("DGEMM", "ZGEMM", "ZHERK"):
        values = dispatch[kernel]
        print(
            f"ozaki_dispatch kernel={kernel} used={values['used']} "
            f"fallback={values['fallback']} unknown={values['unknown']} "
            f"max_bits={values['max_bits']}"
        )
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    print("ozaki_validation=pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
