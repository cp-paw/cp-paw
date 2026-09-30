#!/usr/bin/env python3
"""Measure ionic force differences only after electronic stationarity at each geometry."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess

from cache_parity import digest
from displace_restart import displace, read_records, record_name
import stationarity


def numbers(raw, count):
    fields = raw.split()
    if len(fields) != count or any(not stationarity.NUMBER.fullmatch(x) for x in fields):
        raise ValueError("Malformed force/energy diagnostic")
    values = [float(x.replace("D", "E").replace("d", "e")) for x in fields]
    if not all(math.isfinite(x) for x in values):
        raise ValueError("Nonfinite force/energy diagnostic")
    return values


def force_records(text, natom):
    """Associate each complete force report with its electronic diagnostic step."""
    result, current, step = [], None, None
    for line in text.splitlines():
        if line.startswith(stationarity.HEADER):
            if current is not None or (step is not None and
                                       (not result or result[-1]["step"] != step)):
                raise ValueError("Missing or incomplete force report")
            step = int(line[len(stationarity.HEADER):].strip())
        elif line == "SKALA TOTAL FORCE DIAGNOSTIC":
            if current is not None or step is None or (result and result[-1]["step"] == step):
                raise ValueError("Duplicate or unscoped force report")
            current = {"step": step, "forces": []}
        elif current is not None and line.startswith("TOTAL ENERGY"):
            if "energy" in current:
                raise ValueError("Duplicate total energy")
            current["energy"] = numbers(line[len("TOTAL ENERGY"):], 1)[0]
        elif current is not None and line.startswith("ATOM"):
            fields = line[len("ATOM"):].split()
            if len(fields) != 4 or fields[0] != str(len(current["forces"]) + 1):
                raise ValueError("Invalid or duplicate atom force index")
            current["forces"].append(numbers(" ".join(fields[1:]), 3))
        elif current is not None and line.startswith("NET FORCE"):
            net = numbers(line[len("NET FORCE"):], 3)
            if "energy" not in current or len(current["forces"]) != natom:
                raise ValueError("Incomplete force report")
            if any(not math.isclose(sum(row[i] for row in current["forces"]), net[i],
                                    rel_tol=1e-11, abs_tol=1e-12) for i in range(3)):
                raise ValueError("Net force does not match atom forces")
            result.append(current)
            current = None
    if current is not None or not result or result[-1]["step"] != step:
        raise ValueError("Missing or incomplete force report")
    return result


def geometry(data):
    records = read_records(data)
    sections = {}
    for name in ("ATOMS", "CELL"):
        indices = [i for i, row in enumerate(records) if record_name(row) == name]
        if len(indices) != 1:
            raise ValueError(f"Missing or duplicate {name} section")
        sections[name] = indices[0]
    atom, cell = sections["ATOMS"], sections["CELL"]
    natom, = struct.unpack("<i", records[atom + 1])
    if natom < 1 or len(records[cell + 1]) != 27 * 8:
        raise ValueError("Invalid restart geometry layout")
    positions = [list(struct.unpack(f"<{3 * natom}d", records[atom + offset]))
                 for offset in (3, 4)]
    cells = list(struct.unpack("<27d", records[cell + 1]))
    if not all(math.isfinite(x) for x in [*cells, *positions[0], *positions[1]]):
        raise ValueError("Nonfinite restart geometry")
    return {"natom": natom, "positions": positions, "cells": cells}


def check_geometry(before, after):
    if before["natom"] != after["natom"]:
        raise ValueError("Atom count changed during electronic relaxation")
    # Only the current cell/positions define the stationary geometry. STATIC
    # restarts can legitimately reset the previous propagation time level.
    for old, new in ((before["cells"][:9], after["cells"][:9]),
                     (before["positions"][0], after["positions"][0])):
        if any(not math.isclose(x, y, rel_tol=0, abs_tol=1e-12) for x, y in zip(old, new)):
            raise ValueError("Geometry moved during electronic relaxation")


def displace_geometry(data, atom, axis, delta, rigid_translation=False):
    """Translate all nuclei, or one nucleus, at both coordinate time levels."""
    if not rigid_translation:
        return displace(data, atom, axis, delta)
    for index in range(1, geometry(data)["natom"] + 1):
        data = displace(data, index, axis, delta)
    return data


def force_component(report, atom, axis, rigid_translation=False):
    forces = report["final"]["forces"]
    if axis not in (1, 2, 3) or (not rigid_translation and not 1 <= atom <= len(forces)):
        raise ValueError("Invalid force component")
    return (math.fsum(row[axis - 1] for row in forces) if rigid_translation
            else forces[atom - 1][axis - 1])


def analyze(text, natom, steps, last, residual, commutator):
    trace = stationarity.records(text)
    stationarity.validate(trace)
    if len(trace) != steps or any(b["step"] != a["step"] + 1 for a, b in zip(trace, trace[1:])):
        raise ValueError("Unexpected or nonconsecutive electronic steps")
    detail = stationarity.band_records(text)
    layout = lambda row: [(b["kpoint"], b["spin"], b["band"], b["occupation"])
                          for b in row["bands"]]
    occupations = layout(detail[0])
    if any(layout(row) != occupations for row in detail):
        raise ValueError("Band layout or occupations changed")
    forces = force_records(text, natom)
    if [row["step"] for row in forces] != [row["step"] for row in trace]:
        raise ValueError("Force reports do not match electronic steps")
    result = {"trace": trace, "force_trace": forces, "occupations": occupations,
              "final": forces[-1], "stationary": False}
    try:
        stationarity.validate(trace, residual=residual, commutator=commutator, last=last)
        result["stationary"] = True
    except ValueError as error:
        result["stationarity_failure"] = str(error)
    result["final_energy_span_hartree"] = (max(row["energy"] for row in forces[-last:])
                                            - min(row["energy"] for row in forces[-last:]))
    return result


def compare(center, minus, plus, step, atom, axis, tolerance=None, *, rigid_translation=False):
    if not math.isfinite(step) or step <= 0 or (tolerance is not None and
                                              (not math.isfinite(tolerance) or tolerance <= 0)):
        raise ValueError("Step and optional tolerance must be positive and finite")
    # JSON reloads convert occupation tuples to lists without changing values.
    occupations = [tuple(band) for band in center["occupations"]]
    for row in (center, minus, plus):
        if not row["stationary"]:
            raise ValueError("A geometry has not reached electronic stationarity")
        if [tuple(band) for band in row["occupations"]] != occupations:
            raise ValueError("Changed band layout or occupations across geometries")
    analytic = force_component(center, atom, axis, rigid_translation)
    numeric = -(plus["final"]["energy"] - minus["final"]["energy"]) / (2 * step)
    error = abs(numeric - analytic)
    return {"step_bohr": step, "analytic_hartree_per_bohr": analytic,
            "finite_difference_hartree_per_bohr": numeric,
            "absolute_error_hartree_per_bohr": error,
            "acceptance_limit_hartree_per_bohr": tolerance,
            "passed": None if tolerance is None else error <= tolerance}


def control(args):
    tolerance = getattr(args, "orthogonality_tolerance", None)
    ortho = "" if tolerance is None else f" ORTHOTOL={tolerance:.16e}"
    stress = " STRESS=T" if getattr(args, "stress", False) else ""
    dual = getattr(args, "density_dual", 2.)
    dual_text = "2" if dual == 2. else f"{dual:.16e}"
    return f"""!CONTROL
 !GENERIC TRACE=F DT={args.dt:.16e} NSTEP={args.block_steps} NWRITE=10 START=F
          RSTRTTYPE='STATIC' AUTOCONV=1000 !END
 !DFT TYPE=10
  !SKALA MODEL='model.fun' DEVICE='{args.device}'
   RADIALPOINTS={args.radial_points} LEBEDEVEXACTNESS={args.lebedev_exactness}
   LEBEDEVORIENTATIONS=1 IMAGESHELLS=1 APPLY=T CHECK=T
  !END
 !END
 !FOURIER EPWPSI={args.cutoff:.16e} CDUAL={dual_text} !END
 !CELL MOVE=F FRIC=0.0 M=1.E30 !END
 !PSIDYN SAFEORTHO=T{ortho}{stress} MPSI={args.mass:.16e} MPSICG2={args.mass_g2:.16e}
         FRIC={args.friction:.16e} !END
!END
!EOB
"""


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def run_leg(name, delta, args, inputs, env):
    work = args.output / name
    work.mkdir()
    initial = displace_geometry(inputs["restart"].read_bytes(), args.atom, args.axis, delta,
                                getattr(args, "rigid_translation", False))
    expected = geometry(initial)
    result = {"displacement_bohr": delta, "stationary": False, "blocks": []}
    restart = None
    try:
        for index in range(args.max_blocks):
            block = work / f"block-{index + 1}"
            block.mkdir()
            (block / "si2.cntl").write_text(control(args))
            (block / "model.fun").symlink_to(inputs["model"])
            shutil.copy2(inputs["structure"], block / "si2.strc")
            if restart is None:
                (block / "si2.rstrt").write_bytes(initial)
            else:
                shutil.copy2(restart, block / "si2.rstrt")
            write_json(block / "inputs.json", {key: digest(block / key)
                       for key in ("si2.cntl", "si2.rstrt", "si2.strc", "model.fun")})
            print(f"Running {name}, block {index + 1}: {block}", flush=True)
            with (block / "stdout.log").open("w") as out, (block / "stderr.log").open("w") as err:
                subprocess.run([str(inputs["executable"]), "si2.cntl"], cwd=block,
                               env=env, stdout=out, stderr=err, check=True, timeout=args.timeout)
            restart = block / "si2.rstrt"
            check_geometry(expected, geometry(restart.read_bytes()))
            data = analyze((block / "si2.prot").read_text(), expected["natom"],
                           args.block_steps, args.last, args.residual_tolerance,
                           args.commutator_tolerance)
            if "occupations" in result and result["occupations"] != data["occupations"]:
                raise ValueError("Band layout or occupations changed across restart blocks")
            write_json(block / "results.json", data)
            result["blocks"].append({"directory": str(block), "stationary": data["stationary"],
                                      "final_residual": data["trace"][-1]})
            result.update({key: value for key, value in data.items()
                           if key not in ("trace", "force_trace", "stationarity_failure")})
            if data["stationary"]:
                break
            result["failure"] = data["stationarity_failure"]
        if result["stationary"]:
            result.pop("failure", None)
    except (OSError, ValueError, IndexError, struct.error, subprocess.SubprocessError) as error:
        result.update(stationary=False, failure=str(error))
    write_json(work / "results.json", result)
    print(f"Finished {name}: stationary={result['stationary']}", flush=True)
    return name, result


def add_relaxation_arguments(parser):
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--steps", type=float, nargs="+", default=(0.001, 0.0003))
    parser.add_argument("--block-steps", type=int, default=40)
    parser.add_argument("--max-blocks", type=int, default=3)
    parser.add_argument("--last", type=int, default=5)
    parser.add_argument("--residual-tolerance", type=float, default=1e-6)
    parser.add_argument("--commutator-tolerance", type=float, default=1e-6)
    parser.add_argument("--absolute-tolerance", type=float,
                        help="Omit to measure without claiming passed derivative accuracy")
    parser.add_argument("--radial-points", type=int, default=96)
    parser.add_argument("--lebedev-exactness", type=int, default=17)
    parser.add_argument("--cutoff", type=float, default=20)
    parser.add_argument("--density-dual", type=float, default=2,
                        help="Density/wavefunction cutoff ratio, independently refining the FFT grid")
    parser.add_argument("--dt", type=float, default=5)
    parser.add_argument("--mass", type=float, default=25)
    parser.add_argument("--mass-g2", type=float, default=0.3166286988823056)
    parser.add_argument("--friction", type=float, default=0.05)
    parser.add_argument("--orthogonality-tolerance", type=float,
                        help="Optional tighter PAW constraint solve, from 1e-14 to 1e-8")
    parser.add_argument("--jobs", type=int, default=1, help="Independent serial processes, not MPI ranks")
    parser.add_argument("--timeout", type=float, default=3600, help="Seconds per relaxation block")


def validate_relaxation_arguments(parser, args):
    positive = [args.residual_tolerance, args.commutator_tolerance, args.cutoff, args.density_dual,
                args.dt, args.mass, args.mass_g2, args.friction, args.timeout, *args.steps]
    if args.absolute_tolerance is not None:
        positive.append(args.absolute_tolerance)
    if args.orthogonality_tolerance is not None and not 1e-14 <= args.orthogonality_tolerance <= 1e-8:
        parser.error("Orthogonality tolerance must lie between 1e-14 and 1e-8")
    if (any(not math.isfinite(x) or x <= 0 for x in positive)
            or len(set(args.steps)) != len(args.steps)
            or min(args.jobs, args.block_steps, args.max_blocks, args.last,
                   args.radial_points, args.lebedev_exactness) < 1 or args.last > args.block_steps):
        parser.error("Invalid indices, sizes, relaxation settings or finite-difference parameters")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_relaxation_arguments(parser)
    parser.add_argument("--atom", type=int, default=2)
    parser.add_argument("--axis", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument("--rigid-translation", action="store_true",
                        help="Displace every atom along --axis and test the total force; --atom is unused")
    parser.add_argument("--center-displacement", type=float, default=0)
    args = parser.parse_args()
    validate_relaxation_arguments(parser, args)
    if args.atom < 1 or not math.isfinite(args.center_displacement):
        parser.error("Need a positive atom index and finite displacement")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    geometry(displace_geometry(inputs["restart"].read_bytes(), args.atom, args.axis,
                               args.center_displacement, args.rigid_translation))
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_SCF_DETAIL="1",
               CPPAW_SKALA_DETERMINISTIC="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    provenance = {key: {"path": str(path), "sha256": digest(path)} for key, path in inputs.items()}
    provenance["driver_sha256"] = digest(Path(__file__))
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                               for key, value in vars(args).items()}
    provenance["environment"] = {key: value for key, value in env.items()
                                  if key.startswith("CPPAW_") or key in
                                  ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "CUDA_VISIBLE_DEVICES")}
    write_json(args.output / "provenance.json", provenance)
    tasks = [("center", args.center_displacement), ("center-repeat", args.center_displacement)]
    for index, step in enumerate(args.steps, 1):
        tasks.extend([(f"step-{index}-minus", args.center_displacement - step),
                      (f"step-{index}-plus", args.center_displacement + step)])
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_leg, name, delta, args, inputs, env) for name, delta in tasks]
        legs = dict(future.result() for future in futures)
    summary = {"scope": ("stationary rigid-translation total-force differences at the specified grid"
                         if args.rigid_translation else
                         "stationary fixed-occupation ionic force differences at the specified grid"),
               "quadrature_convergence_certified": False, "legs": legs, "comparisons": [],
               "passed": False, "stationary": all(row["stationary"] for row in legs.values())}
    if summary["stationary"]:
        try:
            if legs["center-repeat"]["occupations"] != legs["center"]["occupations"]:
                raise ValueError("Changed band layout or occupations in repeated center")
            summary["comparisons"] = [compare(legs["center"], legs[f"step-{i}-minus"],
                legs[f"step-{i}-plus"], h, args.atom, args.axis, args.absolute_tolerance,
                rigid_translation=args.rigid_translation)
                for i, h in enumerate(args.steps, 1)]
            summary["center_repeat_energy_difference_hartree"] = abs(
                legs["center"]["final"]["energy"] - legs["center-repeat"]["final"]["energy"])
            summary["center_repeat_force_difference_hartree_per_bohr"] = abs(
                force_component(legs["center"], args.atom, args.axis, args.rigid_translation)
                - force_component(legs["center-repeat"], args.atom, args.axis, args.rigid_translation))
            repeat_consistent = (args.absolute_tolerance is None or (
                summary["center_repeat_force_difference_hartree_per_bohr"] <= args.absolute_tolerance
                and summary["center_repeat_energy_difference_hartree"] / (2 * min(args.steps))
                    <= args.absolute_tolerance))
            summary["passed"] = (None if args.absolute_tolerance is None else
                                 repeat_consistent and all(row["passed"] for row in summary["comparisons"]))
        except ValueError as error:
            summary["failure"] = str(error)
    write_json(args.output / "results.json", summary)
    print(json.dumps({key: value for key, value in summary.items() if key != "legs"}, indent=2))
    return 1 if summary["passed"] is False else 0


if __name__ == "__main__":
    raise SystemExit(main())
