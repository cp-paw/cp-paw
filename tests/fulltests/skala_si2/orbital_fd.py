#!/usr/bin/env python3
"""Measure fixed-occupation orbital derivatives, not stationary ionic forces."""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

from cache_parity import digest
import stationarity


PREFIX = "SKALA ORBITAL ROTATION "
WIDTHS = {"INDICES": 4, "ANGLE": 1, "OCCUPATIONS": 2, "HAMILTONIAN": 2,
          "DERIVATIVE": 1, "PACKED": 1}


def diagnostics(text, direction="subspace"):
    if direction not in ("subspace", "kinetic"):
        raise ValueError("Unknown orbital direction")
    tangent = direction == "kinetic"
    prefix = "SKALA ORBITAL TANGENT " if tangent else PREFIX
    widths = {**WIDTHS, "INDICES": 3, "OCCUPATIONS": 1, "METRIC": 3} if tangent else WIDTHS
    steps = stationarity.records(text)
    stationarity.validate(steps)
    if len(steps) != 1:
        raise ValueError("Expected exactly one electronic step")
    result = {}
    force_section = False
    for line in text.splitlines():
        if line == "SKALA TOTAL FORCE DIAGNOSTIC":
            force_section = True
        elif line.startswith("NET FORCE"):
            force_section = False
        if line.startswith(prefix):
            key, _, raw = line[len(prefix):].partition(" ")
        elif force_section and line.startswith("TOTAL ENERGY "):
            key, raw = "ENERGY", line[len("TOTAL ENERGY "):]
        else:
            continue
        if key in result:
            raise ValueError(f"Duplicate {key}")
        if key == "MODE":
            result[key] = raw.strip()
            if result[key] not in ("REAL", "IMAG"):
                raise ValueError("Invalid rotation mode")
            continue
        fields = raw.split()
        if key not in (*widths, "ENERGY") or len(fields) != widths.get(key, 1):
            raise ValueError(f"Wrong diagnostic layout for {key}")
        if any(not stationarity.NUMBER.fullmatch(x) for x in fields):
            raise ValueError(f"Malformed {key}")
        values = [float(x.replace("D", "E").replace("d", "e")) for x in fields]
        if not all(math.isfinite(x) for x in values):
            raise ValueError(f"Nonfinite {key}")
        result[key] = values
    if set(result) != {*widths, "ENERGY", "MODE"}:
        raise ValueError("Incomplete rotation/energy diagnostics")
    if any(x < 1 or not x.is_integer() for x in result["INDICES"]):
        raise ValueError("Invalid rotation indices")
    if not tangent and result["INDICES"][2] == result["INDICES"][3]:
        raise ValueError("Rotation needs two distinct bands")
    if result["PACKED"] not in ([0.0], [1.0]):
        raise ValueError("Invalid packed-orbital flag")
    if result["PACKED"] == [1.0] and result["MODE"] == "IMAG":
        raise ValueError("Imaginary rotation cannot use packed real orbitals")
    occupations = result["OCCUPATIONS"]
    if min(occupations) < 0 or (occupations[0] == 0 if tangent else
                              occupations[0] == occupations[1]):
        raise ValueError("Probe needs an occupied tangent band or unequal rotation occupations")
    real, imag = result["HAMILTONIAN"]
    if tangent:
        norm, norm_error, orthogonality = result["METRIC"]
        if norm <= 1e-10 or min(norm_error, orthogonality) < 0 or max(norm_error, orthogonality) > 1e-10:
            raise ValueError("Invalid external tangent PAW metric")
        expected = 2 * occupations[0] * real
    else:
        expected = 2 * (occupations[0] - occupations[1])
        expected *= real if result["MODE"] == "REAL" else -imag
    if not math.isclose(result["DERIVATIVE"][0], expected, rel_tol=1e-12, abs_tol=1e-14):
        raise ValueError("Reported derivative disagrees with weighted Hamiltonian")
    result["stationarity"] = steps[0]
    result["DIRECTION"] = direction
    return result


def compare(center, minus, plus, step, absolute_tolerance=None, relative_tolerance=0.0):
    if not math.isfinite(step) or step <= 0:
        raise ValueError("Step must be finite and positive")
    if not math.isfinite(relative_tolerance) or relative_tolerance < 0:
        raise ValueError("Relative tolerance must be finite and nonnegative")
    if absolute_tolerance is not None and (not math.isfinite(absolute_tolerance)
                                          or absolute_tolerance <= 0):
        raise ValueError("Absolute tolerance must be finite and positive")
    for data, offset in ((minus, -step), (plus, step)):
        for key in ("INDICES", "MODE", "OCCUPATIONS", "PACKED", "DIRECTION"):
            if data[key] != center[key]:
                raise ValueError(f"Changed {key} across the finite difference")
        if not math.isclose(data["ANGLE"][0], center["ANGLE"][0] + offset,
                            rel_tol=0, abs_tol=1e-13):
            raise ValueError("Wrong rotation angle")
        if center["DIRECTION"] == "kinetic" and not math.isclose(
                data["METRIC"][0], center["METRIC"][0], rel_tol=1e-12, abs_tol=1e-14):
            raise ValueError("Changed external direction norm across the finite difference")
    analytic = center["DERIVATIVE"][0]
    numeric = (plus["ENERGY"][0] - minus["ENERGY"][0]) / (2 * step)
    error = abs(numeric - analytic)
    limit = None if absolute_tolerance is None else (
        absolute_tolerance + relative_tolerance * abs(analytic))
    return {"step_radians": step, "analytic_hartree": analytic,
            "finite_difference_hartree": numeric, "absolute_error_hartree": error,
            "acceptance_limit_hartree": limit,
            "passed": None if limit is None else error <= limit}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--kpoint", type=int, default=1)
    parser.add_argument("--spin", type=int, default=1)
    parser.add_argument("--direction", choices=("subspace", "kinetic"), default="subspace")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--bands", type=int, nargs=2)
    selection.add_argument("--band", type=int, help="Occupied band for the kinetic tangent")
    parser.add_argument("--mode", choices=("REAL", "IMAG"), default="REAL")
    parser.add_argument("--center-angle", type=float, default=0.02)
    parser.add_argument("--steps", type=float, nargs="+", default=(0.01, 0.003, 0.001))
    parser.add_argument("--radial-points", type=int, default=96)
    parser.add_argument("--lebedev-exactness", type=int, default=17)
    parser.add_argument("--reference-pbe", action="store_true")
    parser.add_argument("--absolute-tolerance", type=float,
                        help="Omit for measurements only, without a pass/fail claim")
    parser.add_argument("--relative-tolerance", type=float, default=0.0)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    if ((args.direction == "kinetic" and args.bands is not None)
            or (args.direction == "subspace" and args.band is not None)):
        parser.error("Use --bands for subspace rotations and --band for kinetic tangents")
    args.bands = (4, 5) if args.bands is None else args.bands
    args.band = 4 if args.band is None else args.band
    if min(args.mpi_ranks, args.kpoint, args.spin, args.band, *args.bands,
           args.radial_points, args.lebedev_exactness) < 1 or args.bands[0] == args.bands[1]:
        parser.error("Indices/grid sizes/ranks must be positive and bands distinct")
    if (not math.isfinite(args.center_angle) or not math.isfinite(args.timeout)
            or args.timeout <= 0 or any(not math.isfinite(h) or h <= 0 for h in args.steps)
            or len(set(args.steps)) != len(args.steps)):
        parser.error("Invalid angle, timeout, or duplicate/nonpositive steps")
    if (not math.isfinite(args.relative_tolerance) or args.relative_tolerance < 0
            or (args.absolute_tolerance is not None and
                (not math.isfinite(args.absolute_tolerance) or args.absolute_tolerance <= 0))):
        parser.error("Invalid tolerances")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode)
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    provenance = {key: {"path": str(path), "sha256": digest(path)}
                  for key, path in inputs.items()}
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                               for key, value in vars(args).items()}
    provenance["environment"] = {key: value for key, value in env.items()
                                  if key.startswith("CPPAW_") or key in
                                  ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    template = Path(__file__).with_name("skala_si2.cntl").read_text()
    replacements = {"START=T": "START=F", "!CELL MOVE=T": "!CELL MOVE=F",
                    "DEVICE='AUTO'": f"DEVICE='{args.device}'",
                    "RADIALPOINTS=200": f"RADIALPOINTS={args.radial_points}",
                    "LEBEDEVEXACTNESS=53": f"LEBEDEVEXACTNESS={args.lebedev_exactness}",
                    "NAME='../si2/stp.cntl'": "NAME='stp.cntl'"}
    if args.reference_pbe:
        replacements["APPLY=T"] = "APPLY=F"
    for old, new in replacements.items():
        if template.count(old) != 1:
            raise ValueError(f"Control template changed: {old}")
        template = template.replace(old, new)
    results = {}

    def run(name, angle):
        work = (args.output / name).resolve()
        work.mkdir()
        (work / "si2.cntl").write_text(template)
        (work / "model.fun").symlink_to(inputs["model"])
        shutil.copy2(inputs["restart"], work / "si2.rstrt")
        shutil.copy2(inputs["structure"], work / "si2.strc")
        shutil.copy2(Path(__file__).parent / "../si2/stp.cntl", work / "stp.cntl")
        indices = [args.kpoint, args.spin, *args.bands]
        variable = "CPPAW_SKALA_ORBITAL_ROTATION"
        if args.direction == "kinetic":
            indices = [args.kpoint, args.spin, args.band]
            variable = "CPPAW_SKALA_ORBITAL_TANGENT"
        env[variable] = " ".join(map(str, indices)) + f" {angle:.17g} {args.mode}"
        command = [str(inputs["executable"]), "si2.cntl"]
        if args.mpi_ranks > 1:
            command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
        print(f"Running {name}, angle {angle:.9g}: {work}", flush=True)
        with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
            subprocess.run(command, cwd=work, env=env, stdout=out, stderr=err,
                           check=True, timeout=args.timeout)
        data = diagnostics((work / "si2.prot").read_text(), args.direction)
        if (data["INDICES"] != indices
                or data["MODE"] != args.mode
                or not math.isclose(data["ANGLE"][0], angle, rel_tol=0, abs_tol=1e-13)):
            raise ValueError("Executable did not apply the requested rotation")
        results[name] = data
        (args.output / "measurements.json").write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n")
        return data

    center = run("center", args.center_angle)
    repeat = run("center-repeat", args.center_angle)
    comparisons = []
    for index, step in enumerate(args.steps):
        minus = run(f"step-{index + 1}-minus", args.center_angle - step)
        plus = run(f"step-{index + 1}-plus", args.center_angle + step)
        comparisons.append(compare(center, minus, plus, step, args.absolute_tolerance,
                                   args.relative_tolerance))
    summary = {"scope": "fixed-occupation orbital derivative, not stationary forces",
               "direction": args.direction,
               "functional": "PBE reference" if args.reference_pbe else "Skala",
               "center_repeat_energy_difference_hartree":
                   abs(center["ENERGY"][0] - repeat["ENERGY"][0]),
               "center_repeat_derivative_difference_hartree":
                   abs(center["DERIVATIVE"][0] - repeat["DERIVATIVE"][0]),
               "comparisons": comparisons,
               "passed": None if args.absolute_tolerance is None else
                   all(row["passed"] for row in comparisons)}
    (args.output / "results.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2, allow_nan=False))
    if summary["passed"] is False:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
