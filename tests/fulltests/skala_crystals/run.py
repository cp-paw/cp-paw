#!/usr/bin/env python3
"""PBE warm start and Skala integration checks for three periodic crystals."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
SPECIES = {"H": (1, "1 1", 1.008), "C": (4, "2 2 1", 12.011),
           "N": (5, "2 2 1", 14.007), "O": (6, "2 2 1", 15.999)}


def read_structure(path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    fields = dict(field.split("=", 1) for field in shlex.split(lines[1]))
    lattice = [float(x) for x in fields["Lattice"].split()]
    atoms = [(row.split()[0], [float(x) for x in row.split()[1:]])
             for row in lines[2:] if row.strip()]
    if len(lattice) != 9 or len(atoms) != count:
        raise ValueError(f"Invalid cell or atom count in {path}")
    for element, xyz in atoms:
        if element not in SPECIES or len(xyz) != 3:
            raise ValueError(f"Invalid atom in {path}: {element} {xyz}")
    if not all(math.isfinite(x) for x in lattice + [x for _, xyz in atoms for x in xyz]):
        raise ValueError(f"Nonfinite geometry in {path}")
    return lattice, atoms


def structure_text(lattice, atoms, kmesh):
    text = ["!STRUCTURE", " !GENERIC LUNIT=1.8897261254578281 STPVERSION='2.0.0' !END",
            f" !KPOINTS DIV={kmesh} {kmesh} {kmesh} !END",
            " !OCCUPATIONS EMPTY=4 NSPIN=1 !END",
            " !LATTICE T=" + " ".join(f"{x:.12e}" for x in lattice) + " !END"]
    for element in dict.fromkeys(element for element, _ in atoms):
        zv, npro, mass = SPECIES[element]
        text.append(f" !SPECIES NAME='{element}_' ZV={zv}. M={mass} NPRO={npro} "
                    f"LRHOX=4 ID='{element.upper()}_.75_6.0' !END")
    for number, (element, xyz) in enumerate(atoms, 1):
        text.append(f" !ATOM NAME='{element}_{number}' R="
                    + " ".join(f"{x:.12e}" for x in xyz) + " !END")
    return "\n".join(text + ["!END", "!EOB", ""])


def control(args, skala):
    steps = args.skala_steps if skala else args.pbe_steps
    functional = " !DFT TYPE=10\n"
    if skala:
        functional += (f"  !SKALA MODEL='model.fun' DEVICE='{args.device}'\n"
                       f"   RADIALPOINTS={args.radial} LEBEDEVEXACTNESS={args.angular}\n"
                       f"   LEBEDEVORIENTATIONS={args.orientations} IMAGESHELLS={args.image_shells}\n"
                       "   APPLY=T CHECK=T !END\n")
    functional += " !END\n"
    # A single-step probe preserves the PBE starting orbitals; longer runs relax them.
    dt = 0.000001 if skala and steps == 1 else args.dt
    return ("!CONTROL\n !FILES\n  !FILE ID='PARMS_STP' NAME='stp.cntl' !END\n !END\n"
            f" !GENERIC TRACE=F DT={dt:.12e} NSTEP={steps} NWRITE={steps} "
            f"START={'F' if skala else 'T'} RSTRTTYPE='STATIC' AUTOCONV=1000 !END\n"
            + functional + f" !FOURIER EPWPSI={args.cutoff} CDUAL=2 !END\n"
            " !PSIDYN MPSI=100 FRIC=0.005\n"
            "  !AUTO FRIC(-)=0.3 FACT(-)=0.97 FRIC(+)=0.3 FACT(+)=1. !END\n"
            " !END\n!END\n!EOB\n")


def last_value(text, label):
    matches = re.findall(r"^" + re.escape(label) + r"\s+([-+0-9.eEdD]+)\s*$", text, re.M)
    if not matches:
        raise ValueError(f"Missing diagnostic: {label}")
    value = float(matches[-1].replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite diagnostic: {label}")
    return value


def execute(executable, work, name):
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    with (work / f"{name}.out").open("w") as out, (work / f"{name}.err").open("w") as err:
        subprocess.run([str(executable), f"{name}.cntl"], cwd=work, env=env,
                       stdout=out, stderr=err, check=True)
    text = (work / f"{name}.prot").read_text()
    if "PROGRAM FINISHED" not in text:
        raise RuntimeError(f"Incomplete run: {work / name}")
    return text


def check_electron_counts(record):
    # The legacy iterative orthogonalizer uses max|S-I| < 1e-8.
    # This bound scales with valence occupations, not frozen-core electrons.
    tolerance = 1.e-8 * max(1., abs(record["OCCUPATION ELECTRONS"]))
    if abs(record["TRACE MINUS OCCUPATIONS"]) > tolerance:
        raise ValueError("PAW overlap trace disagrees with occupations")
    if abs(record["PS GRID MINUS TRACE"]) > 1.e-8:
        raise ValueError("Native pseudo-density grid disagrees with reciprocal norm")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cases", nargs="+", choices=["CO2", "NH3", "urea"],
                        default=["CO2", "NH3", "urea"])
    parser.add_argument("--pbe-steps", type=int, default=180)
    parser.add_argument("--skala-steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=5.0)
    parser.add_argument("--cutoff", type=float, default=40.0)
    parser.add_argument("--kmesh", type=int, default=1)
    parser.add_argument("--radial", type=int, default=200)
    parser.add_argument("--angular", type=int, default=53)
    parser.add_argument("--orientations", type=int, default=1)
    parser.add_argument("--image-shells", type=int, default=1)
    parser.add_argument("--device", choices=["AUTO", "CPU", "CUDA"], default="AUTO")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if min(args.pbe_steps, args.skala_steps, args.radial, args.angular,
           args.kmesh, args.image_shells, args.orientations) < 1:
        parser.error("Step counts and grid settings must be positive")
    if args.angular > 65:
        parser.error("Maximum supported Lebedev exactness is 65")
    if not (math.isfinite(args.dt) and args.dt > 0 and math.isfinite(args.cutoff) and args.cutoff > 0):
        parser.error("Time step and cutoff must be positive and finite")
    executable, model = args.executable.resolve(strict=True), args.model.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=False)
    results = []
    for case in args.cases:
        source = HERE / "structures" / f"{case}-solid.xyz"
        lattice, atoms = read_structure(source)
        work = args.output / case
        work.mkdir()
        shutil.copy2(source, work / source.name)
        shutil.copy2(HERE.parent / "si2/stp.cntl", work / "stp.cntl")
        (work / "model.fun").symlink_to(model)
        for name, skala in [("pbe", False), ("skala", True)]:
            (work / f"{name}.strc").write_text(structure_text(lattice, atoms, args.kmesh))
            (work / f"{name}.cntl").write_text(control(args, skala))
        metadata = {"case": case, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                    "executable": str(executable), "model": str(model),
                    "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
                    "model_sha256": hashlib.sha256(model.read_bytes()).hexdigest(),
                    "settings": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    "scope": "integration smoke; no SCF/grid convergence or cohesive-energy claim"}
        (work / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if args.prepare_only:
            continue
        print(f"{case}: PBE warm start ({args.pbe_steps} steps)", flush=True)
        execute(executable, work, "pbe")
        shutil.copy2(work / "pbe.rstrt", work / "skala.rstrt")
        print(f"{case}: Skala ({args.skala_steps} steps)", flush=True)
        text = execute(executable, work, "skala")
        record = {"case": case}
        for label in ["MODEL XC ENERGY", "COMPOSITE ELECTRONS", "COMPOSITE POSITIVE TAU",
                      "PARTITION VOLUME", "EXACT CELL VOLUME", "PARTITION VOLUME RELATIVE ERROR",
                      "PS RECIPROCAL NORM TRACE", "PS NATIVE GRID ELECTRONS", "PS GRID MINUS TRACE",
                      "PAW VALENCE TRACE", "PAW ALL ELECTRON TRACE", "TRACE MINUS OCCUPATIONS",
                      "OCCUPATION ELECTRONS", "COMPOSITE MINUS TRACE",
                      "SKALA OCCUPIED RESIDUAL RMS", "SKALA OCCUPIED RESIDUAL MAX",
                      "SKALA OCCUPATION COMMUTATOR MAX", "SKALA SCF OVERLAP ERROR",
                      "SKALA HAMILTONIAN HERMITICITY",
                      "GRAD ADJOINT DIFFERENCE", "TAU OPERATOR DIFFERENCE", "ONE-CENTER MATRIX DIFFERENCE"]:
            record[label] = last_value(text, label)
            if "DIFFERENCE" in label and abs(record[label]) > 1.e-8:
                raise ValueError(f"{case}: failed {label}: {record[label]}")
        check_electron_counts(record)
        record["expected_electrons"] = sum({"H": 1, "C": 6, "N": 7, "O": 8}[el] for el, _ in atoms)
        record["electron_quadrature_error"] = record["COMPOSITE ELECTRONS"] - record["expected_electrons"]
        # Record, but do not silently repair, a quadrature or SCF error.
        results.append(record)
        (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
