#!/usr/bin/env python3
"""Relax native SCAN electrons until explicit PAW stationarity criteria hold."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skala_si2"))
import force_stationary_fd as fd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    fd.add_relaxation_arguments(parser)
    parser.add_argument("--electronic-only", action="store_true",
                        help="Require stationary energies without evaluating ionic forces")
    parser.set_defaults(functional="scan", dt=1.5, mass=25., friction=0.4, density_dual=4.,
                        radial_points=64, lebedev_exactness=17, block_steps=5,
                        electronic_warmup_steps=80, max_blocks=20, last=5,
                        residual_tolerance=1e-7, commutator_tolerance=1e-7,
                        orthogonality_tolerance=1e-12)
    args = parser.parse_args()
    fd.validate_relaxation_arguments(parser, args)
    if args.functional != "scan":
        parser.error("Use the general Skala drivers for Skala relaxation")
    inputs = fd.relaxation_inputs(parser, args)
    args.atom, args.axis = 1, 1
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_SCF_DETAIL="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    provenance = {
        "inputs": {key: dict(path=str(path), sha256=fd.digest(path)) for key, path in inputs.items()},
        "driver_sha256": fd.digest(Path(__file__)),
        "arguments": {key: str(value) if isinstance(value, Path) else value
                      for key, value in vars(args).items()},
        "environment": {key: value for key, value in env.items()
                        if key.startswith("CPPAW_") or key in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")},
    }
    fd.write_json(args.output / "provenance.json", provenance)
    _, result = fd.run_leg("relax", 0., args, inputs, env, energy_only=args.electronic_only)
    for key, path in inputs.items():
        if fd.digest(path) != provenance["inputs"][key]["sha256"]:
            raise ValueError(f"Changed source input during relaxation: {key}")
    if result["blocks"]:
        result["restart"] = str(Path(result["blocks"][-1]["directory"]) / "si2.rstrt")
    result["scope"] = "electronic stationarity at the specified PAW basis and quadrature, not cross-code accuracy"
    result["forces_evaluated"] = not args.electronic_only
    fd.write_json(args.output / "results.json", result)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("trace", "force_trace", "energy_trace", "occupations")}, indent=2))
    return 0 if result["stationary"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
