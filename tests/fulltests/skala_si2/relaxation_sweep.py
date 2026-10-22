#!/usr/bin/env python3
"""Compare electronic relaxation parameters at matched residual tolerances."""

import argparse
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace

from cache_parity import digest
import force_stationary_fd as fd


def rank_settings(rows, energy_tolerance):
    accepted = [row for row in rows if row["result"].get("stationary")
                and not row["result"].get("failure")]
    if not accepted:
        return dict(recommended=None, reason="No candidate reached the requested stationarity")
    energies = [row["result"]["final"]["energy"] for row in accepted]
    if any(not math.isfinite(x) for x in energies) or max(energies)-min(energies) > energy_tolerance:
        return dict(recommended=None, reason="Stationary endpoints have inconsistent energies")
    layout = accepted[0]["result"]["occupations"]
    if any(row["result"]["occupations"] != layout for row in accepted):
        return dict(recommended=None, reason="Stationary endpoints have different occupations")
    best = min(accepted, key=lambda row: (row["validated_steps"], row["elapsed_s"]))
    return dict(recommended=best["name"], energy_spread_hartree=max(energies)-min(energies),
                scope="this restart, grid, model, machine and tolerance only",
                selection="fewest validated steps, then elapsed time",
                ground_state_minimum_certified=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    fd.add_relaxation_arguments(parser)
    parser.add_argument("--setting", type=float, nargs=4, action="append", required=True,
                        metavar=("DT", "MASS", "MASS_G2", "FRICTION"))
    parser.add_argument("--energy-tolerance", type=float, default=1e-7)
    parser.set_defaults(block_steps=5, electronic_warmup_steps=20, max_blocks=3, last=5)
    args = parser.parse_args()
    fd.validate_relaxation_arguments(parser, args)
    if args.jobs != 1:
        parser.error("Relaxation candidates must run sequentially (--jobs 1)")
    if (not math.isfinite(args.energy_tolerance) or args.energy_tolerance <= 0 or
            any(not math.isfinite(x) or x <= 0 for setting in args.setting for x in setting)):
        parser.error("Need finite positive settings and energy tolerance")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    hashes = {key: digest(path) for key, path in inputs.items()}
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_DETERMINISTIC="1",
               CPPAW_SKALA_SCF_DETAIL="1", CPPAW_SKALA_TORCH_THREADS="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    fd.write_json(args.output / "provenance.json", dict(
        inputs={key: dict(path=str(path), sha256=hashes[key]) for key, path in inputs.items()},
        arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        environment={key: value for key, value in env.items() if key.startswith("CPPAW_")},
        scope="controlled relaxation search, not universal production defaults"))
    rows = []
    for index, setting in enumerate(args.setting, 1):
        trial = SimpleNamespace(**vars(args))
        trial.dt, trial.mass, trial.mass_g2, trial.friction = setting
        trial.atom, trial.axis = 1, 1
        name = f"setting-{index:02d}"
        start = time.monotonic()
        _, result = fd.run_leg(name, 0., trial, inputs, env)
        validated_steps = (len(result["blocks"])*args.block_steps
                           + len(result.get("electronic_blocks", []))*args.electronic_warmup_steps)
        rows.append(dict(name=name, dt=trial.dt, mass=trial.mass, mass_g2=trial.mass_g2,
                         friction=trial.friction, dt2_over_mass=trial.dt**2/trial.mass,
                         elapsed_s=time.monotonic()-start, validated_steps=validated_steps, result=result))
        fd.write_json(args.output / "results.json", dict(cases=rows, selection=rank_settings(rows, args.energy_tolerance)))
    if any(digest(path) != hashes[key] for key, path in inputs.items()):
        raise ValueError("An executable or scientific input changed during the study")
    print(json.dumps(rank_settings(rows, args.energy_tolerance), indent=2))
    if any("failure" in row["result"] and not row["result"].get("blocks") for row in rows):
        raise SystemExit("One or more settings failed before a valid relaxation block")


if __name__ == "__main__":
    main()
