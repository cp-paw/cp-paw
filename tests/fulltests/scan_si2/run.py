#!/usr/bin/env python3
"""Run native SCAN at a supplied PAW state and audit the active XC operators."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "skala_si2"))
import force_stationary_fd as fd
from force_mode_parity import energy_records
import stationarity
from stress_mode_parity import stress_records


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def control(args, setups):
    flags = "FORCE=F STRESS=F" if args.electronic_only else f"FORCE=T STRESS={'T' if args.stress else 'F'}"
    if "'" in str(setups):
        raise ValueError("The setup path cannot contain a single quote")
    return f"""!CONTROL
 !FILES
  !FILE ID='PARMS_STP' NAME='{setups}' !END
 !END
 !GENERIC TRACE=F DT={args.dt:.16e} NSTEP={args.steps} NWRITE=1 START=F
          RSTRTTYPE='STATIC' AUTOCONV=1000 !END
 !DFT TYPE=10
  !SCAN RADIALPOINTS={args.radial_points} LEBEDEVEXACTNESS={args.lebedev_exactness}
        LEBEDEVORIENTATIONS=1 IMAGESHELLS=1 CHECK=T
  !END
 !END
 !FOURIER EPWPSI={args.cutoff:.16e} CDUAL={args.density_dual:.16e} !END
 !CELL MOVE=F FRIC=0.0 M=1.E30 !END
 !PSIDYN SAFEORTHO=T ORTHOTOL=1.E-12 {flags}
         MPSI={args.mass:.16e} MPSICG2={args.mass_g2:.16e} FRIC={args.friction:.16e} !END
!END
!EOB
"""


def scalar_trace(text, label, steps, *, positive=False):
    values = []
    for line in text.splitlines():
        if line.startswith(label):
            value = fd.numbers(line[len(label):], 1)[0]
            if positive and value <= 0:
                raise ValueError(f"Inactive native SCAN contribution: {label}")
            values.append(value)
    if len(values) != steps:
        raise ValueError(f"Expected {steps} records for {label}, found {len(values)}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("executable", "restart", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    si2 = Path(__file__).resolve().parents[1] / "si2"
    parser.add_argument("--structure", type=Path, default=si2 / "si2.strc")
    parser.add_argument("--setups", type=Path, default=si2 / "stp.cntl")
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1e-6)
    parser.add_argument("--mass", type=float, default=25.)
    parser.add_argument("--mass-g2", type=float, default=0.3166286988823056)
    parser.add_argument("--friction", type=float, default=0.4)
    parser.add_argument("--cutoff", type=float, default=20.)
    parser.add_argument("--density-dual", type=float, default=4.)
    parser.add_argument("--radial-points", type=int, default=64)
    parser.add_argument("--lebedev-exactness", type=int, default=17)
    parser.add_argument("--stress", action="store_true")
    parser.add_argument("--electronic-only", action="store_true")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--residual-tolerance", type=float)
    parser.add_argument("--commutator-tolerance", type=float)
    parser.add_argument("--last", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=1800.)
    args = parser.parse_args()
    if min(args.steps, args.mpi_ranks, args.radial_points, args.lebedev_exactness) < 1:
        parser.error("Step, rank and quadrature counts must be positive")
    for key in ("dt", "mass", "cutoff", "density_dual", "timeout"):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) <= 0:
            parser.error(f"{key} must be positive and finite")
    if not 0 <= args.friction < 1 or not math.isfinite(args.mass_g2) or args.mass_g2 < 0:
        parser.error("Invalid electronic damping or mass preconditioner")
    if args.electronic_only and args.stress:
        parser.error("Electronic-only mode cannot request stress")
    paths = {key: getattr(args, key).resolve(strict=True)
             for key in ("executable", "restart", "structure", "setups")}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    before = fd.geometry(paths["restart"].read_bytes())
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_SCF_DETAIL="1")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    provenance = {
        "inputs": {key: {"path": str(path), "sha256": digest(path)} for key, path in paths.items()},
        "driver_sha256": digest(Path(__file__)),
        "arguments": {key: str(value) if isinstance(value, Path) else value
                      for key, value in vars(args).items()},
    }
    fd.write_json(output / "provenance.json", provenance)
    (output / "si2.cntl").write_text(control(args, paths["setups"]), encoding="ascii")
    shutil.copy2(paths["restart"], output / "si2.rstrt")
    shutil.copy2(paths["structure"], output / "si2.strc")
    command = [str(paths["executable"]), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    with (output / "stdout.log").open("w") as out, (output / "stderr.log").open("w") as err:
        subprocess.run(command, cwd=output, env=env, stdout=out, stderr=err,
                       check=True, timeout=args.timeout)
    text = (output / "si2.prot").read_text()
    if "NATIVE SCAN PAW INTEGRATION" not in text or "SCAN COMPLETE PAW ATOM BLOCKS" not in text:
        raise ValueError("The native SCAN PAW backend was not active")
    trace = stationarity.records(text)
    checks = stationarity.validate(trace, residual=args.residual_tolerance,
                                   commutator=args.commutator_tolerance, last=args.last)
    if len(trace) != args.steps:
        raise ValueError("Native SCAN stopped before the requested electronic steps")
    stationarity.band_records(text)
    reports = energy_records(text) if args.electronic_only else fd.force_records(text, before["natom"])
    if [row["step"] for row in reports] != [row["step"] for row in trace]:
        raise ValueError("Native SCAN energy/force reports do not match the electronic steps")
    operators = {label: scalar_trace(text, label, args.steps, positive=True) for label in (
        "SMOOTH SCALAR OPERATOR L2", "SMOOTH TAU OPERATOR L2", "ONE-CENTER ADJOINT L2")}
    xc = scalar_trace(text, "SCAN XC ENERGY", args.steps)
    core = scalar_trace(text, "SCAN FROZEN-CORE XC REFERENCE", args.steps)
    reconstructed = scalar_trace(text, "SCAN RECONSTRUCTED XC ENERGY", args.steps)
    limits = {label: scalar_trace(text, label, args.steps) for label in (
        "SCAN SPIN ROWS", "SCAN LOW-DENSITY SPIN ROWS", "SCAN NEGATIVE TAU SPIN ROWS",
        "SCAN PAULI-BOUND VIOLATING SPIN ROWS", "SCAN PAULI-BOUND VIOLATING ELECTRONS",
        "SCAN ABSOLUTE QUADRATURE ELECTRONS", "SCAN INTEGRATED PAULI DEFICIT",
        "SCAN INTEGRATED NEGATIVE TAU DEFICIT")}
    if any(not math.isclose(raw-ref, net, rel_tol=1e-12, abs_tol=1e-12)
           for raw, ref, net in zip(reconstructed, core, xc)):
        raise ValueError("Native SCAN frozen-core energy decomposition is inconsistent")
    stress = stress_records(text) if args.stress else None
    fd.check_geometry(before, fd.geometry((output / "si2.rstrt").read_bytes()))
    for key, path in paths.items():
        if digest(path) != provenance["inputs"][key]["sha256"]:
            raise ValueError(f"Changed source input during verification: {key}")
    result = {"passed": True, "scope": checks["scope"], "stationarity": checks,
              "native_scan_xc_hartree": xc, "reconstructed_xc_hartree": reconstructed,
              "frozen_core_xc_hartree": core, "operators": operators, "final": reports[-1],
              "stress": stress, "field_limits": limits}
    fd.write_json(output / "results.json", result)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
