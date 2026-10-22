#!/usr/bin/env python3
"""Check complete-block batching before changing native Skala inference."""

import argparse
import json
import math
import os
from pathlib import Path

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import torch
from model_precision import features, evaluate


def blocks(count, device):
    original = features(device)
    result = []
    for atom in range(count):
        source = atom % 2
        first, last = source*8, source*8 + 8-atom%3
        delta = torch.tensor([2.3*atom, -.2*atom, .4*atom], dtype=torch.float64, device=device)
        current = {}
        for key, value in original.items():
            if key in ("density", "grad", "kin"):
                current[key] = value[..., first:last].clone()
            elif key in ("grid_coords", "grid_weights", "atomic_grid_weights"):
                current[key] = value[first:last].clone()
            elif key == "coarse_0_atomic_coords":
                current[key] = value[source:source+1].clone()
            elif key == "atomic_grid_sizes":
                current[key] = torch.tensor([last-first], dtype=torch.int64, device=device)
            else:
                current[key] = value[:last-first].clone()
        current["grid_coords"] += delta
        current["coarse_0_atomic_coords"] += delta
        result.append(current)
    return result


def combine(parts):
    return {key: (parts[0][key].new_empty((max(len(p["grid_weights"]) for p in parts), 0))
                  if key == "atomic_grid_size_bound_shape" else
                  torch.cat([p[key] for p in parts], dim=-1 if key in ("density", "grad", "kin") else 0))
            for key in parts[0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--blocks", type=int, nargs="+", default=[2, 4, 8])
    parser.add_argument("--tolerance", type=float, default=2e-7)
    args = parser.parse_args()
    if args.output.exists() or min(args.blocks) < 1 or not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("Need a new output file, positive block counts and finite tolerance")
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    model = torch.jit.load(str(args.model), map_location=args.device).eval().to(args.device)
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    rows = []
    for count in args.blocks:
        parts = blocks(count, args.device)
        separate = [evaluate(model, p, derivatives=True) for p in parts]
        energy, derivatives = evaluate(model, combine(parts), derivatives=True)
        differences = {"energy": abs(energy-math.fsum(row[0] for row in separate))}
        if any(set(row[1]) != set(derivatives) for row in separate):
            raise ValueError("Batching changed derivative feature set")
        for key, value in derivatives.items():
            expected = torch.cat([row[1][key] for row in separate],
                                 dim=-1 if key in ("density", "grad", "kin") else 0)
            differences[key] = (value-expected).abs().max().item()
        rows.append(dict(blocks=count, differences=differences,
                         passed=all(math.isfinite(x) and x <= args.tolerance for x in differences.values())))
    from cache_parity import digest
    result = dict(model=str(args.model.resolve()), sha256=digest(args.model), torch=torch.__version__,
                  device=args.device, tolerance=args.tolerance, cases=rows,
                  scope="synthetic complete-block separability, not physical reference validation",
                  passed=all(row["passed"] for row in rows))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit("FAIL: batching is not established as a numerically equivalent execution mode")


if __name__ == "__main__":
    main()
