#!/usr/bin/env python3
"""Validate CHECK=T electronic diagnostics, optionally requiring stationarity."""

import argparse
import json
import math
from pathlib import Path
import re


HEADER = "SKALA ELECTRONIC STATIONARITY STEP"
LABELS = {
    "SKALA OCCUPIED RESIDUAL RMS": "rms",
    "SKALA OCCUPIED RESIDUAL MAX": "maximum",
    "SKALA OCCUPATION COMMUTATOR MAX": "commutator",
    "SKALA SCF OVERLAP ERROR": "overlap",
    "SKALA HAMILTONIAN HERMITICITY": "hermiticity",
}
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


def records(text):
    if "PROGRAM FINISHED" not in text:
        raise ValueError("Calculation has not finished normally")
    result = []
    current = None
    for line in text.splitlines():
        if line.startswith(HEADER):
            step = line[len(HEADER):].strip()
            if not step.isdigit():
                raise ValueError("Invalid electronic diagnostic step")
            current = {"step": int(step)}
            result.append(current)
        for label, key in LABELS.items():
            if not line.startswith(label):
                continue
            raw = line[len(label):].strip()
            if current is None or key in current or not NUMBER.fullmatch(raw):
                raise ValueError(f"Invalid or duplicate {label}")
            value = float(raw.replace("D", "E").replace("d", "e"))
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Non-finite or negative {label}")
            current[key] = value
    if not result:
        raise ValueError("No electronic stationarity diagnostics; enable CHECK=T")
    for record in result:
        if set(record) != {"step", *LABELS.values()}:
            raise ValueError(f"Incomplete electronic diagnostics at step {record['step']}")
        if record["rms"] > record["maximum"] + 1e-12 * max(1.0, record["maximum"]):
            raise ValueError("Occupied RMS exceeds the maximum residual")
    return result


def validate(data, *, overlap=1e-8, hermiticity=1e-10,
             residual=None, commutator=None, last=1):
    if not data or last < 1 or last > len(data):
        raise ValueError("Not enough diagnostic steps for the requested final window")
    if (residual is None) != (commutator is None):
        raise ValueError("Supply both residual and occupation-commutator tolerances")
    limits = [overlap, hermiticity]
    if residual is not None:
        limits.extend([residual, commutator])
    if any(not math.isfinite(x) or x <= 0 for x in limits):
        raise ValueError("Tolerances must be positive and finite")
    # Metric/operator failures at earlier steps must not be hidden by the last one.
    for record in data:
        for key, limit in (("overlap", overlap), ("hermiticity", hermiticity)):
            if record[key] > limit:
                raise ValueError(f"Step {record['step']}: {key}={record[key]:.6e} > {limit:.6e}")
    if residual is not None:
        for record in data[-last:]:
            for key, limit in (("maximum", residual), ("commutator", commutator)):
                if record[key] > limit:
                    raise ValueError(f"Step {record['step']}: {key}={record[key]:.6e} > {limit:.6e}")
    return {"scope": "stationarity" if residual is not None else "diagnostic consistency",
            "steps": len(data), "final": data[-1],
            "ground_state_minimum_certified": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protocol", type=Path)
    parser.add_argument("--residual-tolerance", type=float)
    parser.add_argument("--commutator-tolerance", type=float)
    parser.add_argument("--overlap-tolerance", type=float, default=1e-8)
    parser.add_argument("--hermiticity-tolerance", type=float, default=1e-10)
    parser.add_argument("--last", type=int, default=1)
    args = parser.parse_args()
    try:
        summary = validate(records(args.protocol.read_text()),
                           overlap=args.overlap_tolerance,
                           hermiticity=args.hermiticity_tolerance,
                           residual=args.residual_tolerance,
                           commutator=args.commutator_tolerance, last=args.last)
    except (OSError, ValueError) as error:
        parser.exit(1, f"TEST FAILED: {error}\n")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
