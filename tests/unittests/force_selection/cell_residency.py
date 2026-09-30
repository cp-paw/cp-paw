#!/usr/bin/env python3
"""Check moving-cell host/device ownership with conventional Si2 and a GPU profile build."""

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "fulltests/skala_si2"))
import force_mode_parity as modes
import force_stationary_fd as fd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("Rank count must be positive")
    binary = args.executable.resolve(strict=True)
    binary_sha256 = fd.digest(binary)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixture = Path(__file__).resolve().parents[2] / "fulltests/si2"
    template = (fixture / "si2.cntl").read_text()
    template = template.replace("DT=10.0", "DT=1.0").replace("NSTEP=180", "NSTEP=4")
    template = template.replace("!CONTROL\n", "!CONTROL\n !ANALYSE !TRA FORCE=T E=T !END !END\n")
    env = {key: value for key, value in os.environ.items() if not key.startswith("CPPAW_")}
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               CPPAW_ACCEL_PROFILE="1", CPPAW_ACCEL_PROFILE_FILE="ownership")
    command = [str(binary), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    fd.write_json(root / "provenance.json", dict(binary_sha256=binary_sha256,
        driver_sha256=fd.digest(Path(__file__)), mpi_ranks=args.mpi_ranks,
        scope="moving-cell implementation parity, not physical convergence"))
    results = {}
    for name in ("seed", "off", "transfer", "resident"):
        work = root / name
        work.mkdir()
        for filename in ("si2.strc", "stp.cntl"):
            shutil.copy2(fixture / filename, work / filename)
        control = template
        if name != "seed":
            control = control.replace("START=t", "START=f")
            control = control.replace("!PSIDYN", "!CELL MOVE=T FRIC=0.0 M=1.E6 !END\n !PSIDYN")
            shutil.copy2(root / "seed/si2.rstrt", work / "si2.rstrt")
        (work / "si2.cntl").write_text(control)
        fd.write_json(work / "inputs.json", {p.name: fd.digest(p) for p in work.iterdir()
                                              if p.is_file()})
        run_env = dict(env, CPPAW_GPU_MODE="off" if name == "seed" else name)
        print(f"Running {name}", flush=True)
        with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
            subprocess.run(command, cwd=work, env=run_env, stdout=out, stderr=err,
                           timeout=600, check=True)
        if "PROGRAM FINISHED" not in (work / "si2.prot").read_text():
            raise ValueError(f"{name}: missing normal termination")
        restart = (work / "si2.rstrt").read_bytes()
        results[name] = dict(geometry=fd.geometry(restart), waves=modes.wave_records(restart))
    reference = results["off"]
    movement = max(abs(a-b) for a, b in zip(reference["geometry"]["cells"][:9],
                                            results["seed"]["geometry"]["cells"][:9]))
    if movement < 1e-10:
        raise ValueError("Cell did not move measurably")
    differences = {}
    for name in ("transfer", "resident"):
        current = results[name]
        differences[name] = dict(
            cell=max(abs(a-b) for a, b in zip(reference["geometry"]["cells"],
                                             current["geometry"]["cells"])),
            waves_and_lambda=modes.wave_difference(reference["waves"], current["waves"],
                                                    cell_tolerance=1e-9))
        for label, suffix in (("energy", "e"), ("force", "f")):
            differences[name][label] = modes.trajectory_difference(
                (root / "off" / f"si2_{suffix}.tra").read_bytes(),
                (root / name / f"si2_{suffix}.tra").read_bytes())
        if max(differences[name].values()) > 1e-8:
            raise ValueError(f"{name}: moving-cell mismatch {differences[name]}")
    counts = {}
    profiles = list((root / "resident").glob("ownership*.csv"))
    if len(profiles) != args.mpi_ranks:
        raise ValueError("Missing rank-specific GPU profiles")
    for profile in profiles:
        with profile.open() as stream:
            for row in csv.DictReader(stream):
                counts[row["op"]] = counts.get(row["op"], 0) + int(row["calls"])
    for label in ("ACC_COPY_SETUP_PSIM_IN", "ACC_COPY_PROP_PSIM_HOST_OUT"):
        if counts.get(label, 0) < 1:
            raise ValueError(f"Ownership transition not exercised: {label}")
    if fd.digest(binary) != binary_sha256:
        raise ValueError("Executable changed during test")
    fd.write_json(root / "results.json", dict(passed=True, differences=differences,
        cell_motion=movement, ownership_calls={key: counts[key] for key in
            ("ACC_COPY_SETUP_PSIM_IN", "ACC_COPY_PROP_PSIM_HOST_OUT")}))
    print(differences, flush=True)


if __name__ == "__main__":
    main()
