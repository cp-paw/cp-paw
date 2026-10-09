#!/usr/bin/env python3
"""Measure strain derivatives after independently relaxing every strained cell."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import struct

from deform_restart import deform, write_records
import force_mode_parity as modes
import force_stationary_fd as fd
from stress_mode_parity import stress_records


def direction(name):
    matrix = [[0. for _ in range(3)] for _ in range(3)]
    if name == "isotropic":
        for i in range(3):
            matrix[i][i] = 1.
    elif name in ("xx", "yy", "zz", "xy", "xz", "yz"):
        i, j = ("xyz".index(letter) for letter in name)
        matrix[i][j] = matrix[j][i] = 1. if i == j else 0.5
    else:
        raise ValueError("Unknown strain direction")
    return matrix


def basis_signature(data):
    """Hash integer G vectors, k points and band/spin metadata, not cell or waves."""
    records = modes.wave_records(data)
    modes.wave_difference(records, records)
    nk, ns = struct.unpack_from("<2i", records[0])
    nw, = struct.unpack_from("<i", records[0], 80)
    metadata = [records[0][:8] + records[0][80:]]
    index = 1
    for _ in range(nk):
        _, ng, ndim, nb, real = struct.unpack("<8s4i", records[index])
        metadata.extend(records[index:index+2])
        index += 2 + nw * ns * ((nb + 1) // 2 if real else nb)
        metadata.append(records[index])
        _, _, _, _, nl = struct.unpack("<8s4i", records[index])
        index += 1 + ns * nl
    return hashlib.sha256(write_records(metadata)).hexdigest()


def grid_signature(text):
    labels = ("#(G-VECTORS FOR DENSITY)", "GRID POINTS")
    values = []
    for label in labels:
        rows = [line for line in text.splitlines() if line.startswith(label)]
        if len(rows) != 1:
            raise ValueError("Missing or duplicate grid-size diagnostic")
        value = rows[0].split()[-1]
        if not value.isdigit() or int(value) < 1:
            raise ValueError("Invalid grid-size diagnostic")
        values.append(int(value))
    return values


def compare(center, minus, plus, step, name, tolerance=None):
    if not math.isfinite(step) or step <= 0 or (tolerance is not None and
            (not math.isfinite(tolerance) or tolerance <= 0)):
        raise ValueError("Step and optional tolerance must be positive and finite")
    occupations = [tuple(band) for band in center["occupations"]]
    for row in (center, minus, plus):
        if not row["stationary"]:
            raise ValueError("A strained geometry is not electronically stationary")
        if [tuple(band) for band in row["occupations"]] != occupations:
            raise ValueError("Changed occupations across strained geometries")
        if any(row[key] != center[key] for key in ("basis_signature", "grid_signature")):
            raise ValueError("Changed basis or grid across strained geometries")
    tensor = center["final_stress"]["tensor"]
    strain = direction(name)
    analytic = math.fsum(tensor[i][j] * strain[i][j] for i in range(3) for j in range(3))
    numeric = (plus["final"]["energy"]-minus["final"]["energy"])/(2*step)
    error = abs(numeric-analytic)
    return dict(direction=name, step=step, analytic_hartree=analytic,
                finite_difference_hartree=numeric, absolute_error_hartree=error,
                energy_span_sensitivity_hartree=(minus["final_energy_span_hartree"]
                    + plus["final_energy_span_hartree"])/(2*step),
                acceptance_limit_hartree=tolerance,
                passed=None if tolerance is None else error <= tolerance)


def run_leg(name, components, args, inputs, env):
    seed = args.output / (name + "-initial.rstrt")
    seed.write_bytes(deform(inputs["restart"].read_bytes(), components))
    changed_inputs = dict(inputs, restart=seed)
    _, result = fd.run_leg(name, 0., args, changed_inputs, env)
    result.pop("displacement_bohr")
    result["strain"] = components
    result["initial_restart_sha256"] = fd.digest(seed)
    try:
        basis, grid = basis_signature(inputs["restart"].read_bytes()), None
        electronic_blocks = result.get("electronic_blocks", [])
        for block in electronic_blocks + result["blocks"]:
            work = Path(block["directory"])
            text = (work / "si2.prot").read_text()
            current = grid_signature(text)
            if grid is not None and current != grid:
                raise ValueError("Grid size changed across relaxation blocks")
            grid = current
            if basis_signature((work / "si2.rstrt").read_bytes()) != basis:
                raise ValueError("Plane-wave basis changed during relaxation")
            if block in electronic_blocks:
                continue
            trace = stress_records(text)
            fd.write_json(work / "stress.json", trace)
            result["final_stress"] = trace[-1]
            result["final_stress_span"] = max(
                max(row["tensor"][i][j] for row in trace[-args.last:])
                - min(row["tensor"][i][j] for row in trace[-args.last:])
                for i in range(3) for j in range(3))
        if not result["blocks"]:
            raise ValueError(result.get("failure", "No completed relaxation block"))
        if result["final_stress"]["step"] != result["final"]["step"]:
            raise ValueError("Stress and energy belong to different evaluations")
        result.update(basis_signature=basis, grid_signature=grid)
    except (OSError, ValueError, KeyError, struct.error) as error:
        result.update(stationary=False, failure=str(error))
    fd.write_json(args.output / name / "results.json", result)
    return name, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    fd.add_relaxation_arguments(parser)
    parser.add_argument("--directions", nargs="+", choices=("isotropic", "xx", "yy", "zz", "xy", "xz", "yz"),
                        default=("isotropic", "xx", "xy"))
    parser.set_defaults(steps=(3e-5, 1e-5), residual_tolerance=1e-7,
                        commutator_tolerance=1e-7, orthogonality_tolerance=1e-12)
    args = parser.parse_args()
    fd.validate_relaxation_arguments(parser, args)
    if len(set(args.directions)) != len(args.directions):
        parser.error("Strain directions must be unique")
    args.stress, args.atom, args.axis = True, 1, 1
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    fd.geometry(inputs["restart"].read_bytes())
    basis_signature(inputs["restart"].read_bytes())
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_SCF_DETAIL="1",
               CPPAW_SKALA_DETERMINISTIC="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    provenance = {key: dict(path=str(path), sha256=fd.digest(path)) for key, path in inputs.items()}
    provenance.update(driver_sha256=fd.digest(Path(__file__)),
        arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        environment={key: value for key, value in env.items() if key.startswith("CPPAW_") or key in
                     ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")})
    fd.write_json(args.output / "provenance.json", provenance)
    tasks = [("center", [0.]*9), ("center-repeat", [0.]*9)]
    for name in args.directions:
        matrix = direction(name)
        for index, step in enumerate(args.steps, 1):
            for sign, value in (("minus", -step), ("plus", step)):
                tasks.append((f"{name}-{index}-{sign}", [value*x for row in matrix for x in row]))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_leg, name, components, args, inputs, env) for name, components in tasks]
        legs = dict(future.result() for future in futures)
    for key, path in inputs.items():
        if fd.digest(path) != provenance[key]["sha256"]:
            raise ValueError(f"Changed input during strain verification: {key}")
    summary = dict(scope="stationary strain derivatives at fixed fractional nuclei and finite resolution",
                   quadrature_convergence_certified=False, legs=legs, comparisons=[], passed=False,
                   stationary=all(row["stationary"] for row in legs.values()))
    if summary["stationary"]:
        try:
            compare(legs["center"], legs["center-repeat"], legs["center-repeat"],
                    args.steps[0], args.directions[0])
            summary["comparisons"] = [compare(legs["center"], legs[f"{name}-{i}-minus"],
                legs[f"{name}-{i}-plus"], step, name, args.absolute_tolerance)
                for name in args.directions for i, step in enumerate(args.steps, 1)]
            repeat_energy = abs(legs["center"]["final"]["energy"]-legs["center-repeat"]["final"]["energy"])
            repeat_stress = max(abs(a-b) for first, second in
                zip(legs["center"]["final_stress"]["tensor"], legs["center-repeat"]["final_stress"]["tensor"])
                for a, b in zip(first, second))
            summary.update(center_repeat_energy_difference_hartree=repeat_energy,
                           center_repeat_stress_difference_hartree=repeat_stress)
            summary["passed"] = (None if args.absolute_tolerance is None else
                all(row["passed"] for row in summary["comparisons"])
                and max(repeat_stress, repeat_energy/(2*min(args.steps))) <= args.absolute_tolerance)
        except ValueError as error:
            summary["failure"] = str(error)
    fd.write_json(args.output / "results.json", summary)
    print(json.dumps({key: value for key, value in summary.items() if key != "legs"}, indent=2))
    return 1 if summary["passed"] is False else 0


if __name__ == "__main__":
    raise SystemExit(main())
