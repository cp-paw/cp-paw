#!/usr/bin/env python3
"""Check fixed-geometry force omission without accepting an incomplete force."""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
from types import SimpleNamespace

from cache_parity import digest
from displace_restart import read_records, record_name
import force_stationary_fd as fd
import stationarity


def energy_records(text):
    """Require one complete energy-only report per electronic step."""
    trace = stationarity.records(text)
    expected = [row["step"] for row in trace]
    result, current, step = [], None, None
    for line in text.splitlines():
        if line.startswith(stationarity.HEADER):
            if current is not None or (step is not None and
                                       (not result or result[-1]["step"] != step)):
                raise ValueError("Missing energy-only report")
            step = int(line[len(stationarity.HEADER):].strip())
        elif line == "SKALA TOTAL ENERGY DIAGNOSTIC":
            if step is None or current is not None or (result and result[-1]["step"] == step):
                raise ValueError("Duplicate or unscoped energy-only report")
            current = {"step": step}
        elif current is not None and line.startswith("TOTAL ENERGY"):
            if "energy" in current:
                raise ValueError("Duplicate total energy")
            current["energy"] = fd.numbers(line[len("TOTAL ENERGY"):], 1)[0]
        elif line == "NUCLEAR FORCES NOT CALCULATED":
            if current is None or "energy" not in current:
                raise ValueError("Incomplete energy-only report")
            result.append(current)
            current = None
        elif line == "SKALA TOTAL FORCE DIAGNOSTIC" or line.startswith("SKALA FORCE ATOM"):
            raise ValueError("Unexpected force output in energy-only calculation")
        elif current is not None and (line.startswith("ATOM ") or line.startswith("NET FORCE")):
            raise ValueError("Unexpected forces in energy-only report")
    if current is not None or [row["step"] for row in result] != expected:
        raise ValueError("Energy-only reports do not match electronic steps")
    return result


def wave_records(data):
    records = read_records(data)
    indices = [i for i, row in enumerate(records) if record_name(row) == "WAVES"]
    if len(indices) != 1:
        raise ValueError("Missing or duplicate WAVES restart section")
    start = indices[0]
    count, = struct.unpack_from("<i", records[start])
    if count < 1 or start + count >= len(records):
        raise ValueError("Invalid WAVES record count")
    return records[start + 1:start + count + 1]


def wave_difference(a, b, *, cell_tolerance=0.):
    """Compare numeric wave/Lambda payloads, requiring identical binary metadata."""
    if not math.isfinite(cell_tolerance) or cell_tolerance < 0:
        raise ValueError("Cell tolerance must be finite and nonnegative")
    if (len(a) != len(b) or not a or len(a[0]) != 84 or len(b[0]) != 84
            or a[0][:8] != b[0][:8] or a[0][80:] != b[0][80:]):
        raise ValueError("Changed electronic restart metadata")
    cells = [struct.unpack_from("<9d", header, 8) for header in (a[0], b[0])]
    if not all(math.isfinite(x) for cell in cells for x in cell):
        raise ValueError("Nonfinite electronic restart cell")
    if max(abs(x-y) for x, y in zip(*cells)) > cell_tolerance:
        raise ValueError("Changed electronic restart cell")
    nk, ns = struct.unpack_from("<2i", a[0])
    nw, = struct.unpack_from("<i", a[0], 80)
    if min(nk, ns, nw) < 1:
        raise ValueError("Invalid electronic restart dimensions")
    index, maximum = 1, 0.

    def metadata(size):
        nonlocal index
        if index >= len(a) or len(a[index]) != size or a[index] != b[index]:
            raise ValueError("Changed electronic restart metadata")
        record = a[index]
        index += 1
        return record

    def numeric(count):
        nonlocal index, maximum
        if index >= len(a) or len(a[index]) != 8 * count or len(b[index]) != 8 * count:
            raise ValueError("Invalid electronic restart payload")
        x, y = (struct.unpack(f"<{count}d", records[index]) for records in (a, b))
        if not all(math.isfinite(v) for v in (*x, *y)):
            raise ValueError("Nonfinite electronic restart")
        maximum = max(maximum, max(abs(u - v) for u, v in zip(x, y)))
        index += 1

    for _ in range(nk):
        key, ng, ndim, nb, real = struct.unpack("<8s4i", metadata(24))
        if key.strip() != b"PSI" or min(ng, ndim, nb) < 1 or real not in (0, 1):
            raise ValueError("Invalid PSI restart header")
        metadata(24 + 12 * ng)
        nbh = (nb + 1) // 2 if real else nb
        for _ in range(nw * ns * nbh):
            numeric(2 * ng * ndim)
        key, n, nd, nspin, nl = struct.unpack("<8s4i", metadata(24))
        if key.strip() != b"LAMBDA" or (n, nd, nspin) != (nb, ndim, ns) or nl < 1:
            raise ValueError("Invalid LAMBDA restart header")
        for _ in range(ns * nl):
            numeric(2 * nb * nb)
    if index != len(a):
        raise ValueError("Unparsed electronic restart records")
    return maximum


def scalars(text, count):
    labels = ("MODEL XC ENERGY", "COMPOSITE ELECTRONS", "COMPOSITE POSITIVE TAU",
              "SMOOTH RHO ADJOINT L2", "SMOOTH GRAD ADJOINT L2", "SMOOTH TAU ADJOINT L2",
              "ONE-CENTER ADJOINT L2", "ONE-CENTER MATRIX DIFFERENCE",
              "SMOOTH SCALAR OPERATOR L2", "SMOOTH TAU OPERATOR L2")
    result = {label: [] for label in labels}
    for line in text.splitlines():
        for label in labels:
            if line.startswith(label + " "):
                result[label].append(fd.numbers(line[len(label):], 1)[0])
    if any(len(values) != count for values in result.values()):
        raise ValueError("Incomplete electronic adjoint diagnostics")
    return result


def trajectory_difference(first, second):
    a, b = read_records(first), read_records(second)
    if len(a) != len(b) or not a:
        raise ValueError("Changed or empty trajectory")
    difference = 0.
    for x, y in zip(a, b):
        if len(x) < 16 or x[:16] != y[:16]:
            raise ValueError("Changed trajectory step, time or size")
        _, time, count = struct.unpack("<idi", x[:16])
        if not math.isfinite(time) or count < 1 or len(x) != 16 + 8 * count or len(y) != len(x):
            raise ValueError("Invalid trajectory payload")
        u, v = (struct.unpack(f"<{count}d", row[16:]) for row in (x, y))
        if not all(math.isfinite(value) for value in (*u, *v)):
            raise ValueError("Nonfinite trajectory")
        difference = max(difference, max(abs(p - q) for p, q in zip(u, v)))
    return difference


def compare(on, off, tolerance):
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Tolerance must be positive and finite")
    if len(on["energies"]) != len(off["energies"]) or not on["energies"]:
        raise ValueError("Changed or empty energy trace")
    count = len(on["energies"])
    if set(on["scalars"]) != set(off["scalars"]) or not on["scalars"]:
        raise ValueError("Changed or empty diagnostic layout")
    for data in (on, off):
        if len(data["trace"]) != count or len(data["bands"]) != count:
            raise ValueError("Changed diagnostic step count")
        steps = [row["step"] for row in data["trace"]]
        if steps != [row["step"] for row in data["energies"]] or \
                steps != [row["step"] for row in data["bands"]]:
            raise ValueError("Mismatched diagnostic steps")
    if [row["step"] for row in on["trace"]] != [row["step"] for row in off["trace"]]:
        raise ValueError("Changed electronic steps")
    differences = {}
    for label in on["scalars"]:
        a, b = on["scalars"][label], off["scalars"][label]
        if len(a) != count or len(b) != count:
            raise ValueError("Changed diagnostic layout")
        differences[label] = max(abs(x - y) for x, y in zip(a, b))
    differences["TOTAL ENERGY"] = max(abs(a["energy"] - b["energy"])
                                       for a, b in zip(on["energies"], off["energies"]))
    for key in stationarity.LABELS.values():
        differences[key] = max(abs(a[key] - b[key]) for a, b in zip(on["trace"], off["trace"]))
    band_difference = 0.
    for a, b in zip(on["bands"], off["bands"]):
        if len(a["bands"]) != len(b["bands"]) or a["step"] != b["step"]:
            raise ValueError("Changed band layout")
        for x, y in zip(a["bands"], b["bands"]):
            if any(x[key] != y[key] for key in ("kpoint", "spin", "band", "occupation")):
                raise ValueError("Changed band identity or occupation")
            band_difference = max(band_difference, *(abs(x[key] - y[key])
                                  for key in ("expectation", "residual", "commutator")))
    differences["band_diagnostics"] = band_difference
    for label, value in differences.items():
        if not math.isfinite(value) or value > tolerance:
            raise ValueError(f"{label}: force-mode difference {value} > {tolerance}")
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "model", "restart", "structure", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--device", choices=("CPU", "CUDA"), default="CPU")
    parser.add_argument("--gpu-mode", choices=("off", "transfer", "resident"), default="off")
    parser.add_argument("--source-cache-mib", type=int, default=0)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    if args.source_cache_mib < 0 or args.mpi_ranks < 1:
        parser.error("Need a nonnegative cache budget and positive rank count")
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("Tolerance must be positive and finite")
    inputs = {key: getattr(args, key).resolve(strict=True)
              for key in ("executable", "model", "restart", "structure")}
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, CPPAW_GPU_MODE=args.gpu_mode, CPPAW_SKALA_DETERMINISTIC="1",
               CPPAW_SKALA_SCF_DETAIL="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               CPPAW_SKALA_SOURCE_CACHE_MB=str(args.source_cache_mib),
               CPPAW_SKALA_PARTITION_CACHE_MB="256")
    for key in ("CPPAW_SKALA_ORBITAL_ROTATION", "CPPAW_SKALA_ORBITAL_TANGENT"):
        env.pop(key, None)
    settings = SimpleNamespace(dt=0.001, block_steps=3, device=args.device,
                               radial_points=96, lebedev_exactness=17, cutoff=20.,
                               mass=25., mass_g2=0.3166286988823056, friction=0.4,
                               orthogonality_tolerance=1e-12)
    template = fd.control(settings).replace("!CONTROL\n", "!CONTROL\n !ANALYSE !TRA FORCE=T !END !END\n")
    provenance = {key: {"path": str(path), "sha256": digest(path)} for key, path in inputs.items()}
    provenance["arguments"] = {key: str(value) if isinstance(value, Path) else value
                                for key, value in vars(args).items()}
    provenance["environment"] = {k: v for k, v in env.items() if k.startswith("CPPAW_") or
                                 k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")}
    fd.write_json(args.output / "provenance.json", provenance)
    command = [str(inputs["executable"]), "si2.cntl"]
    if args.mpi_ranks > 1:
        command = [args.mpiexec, "-np", str(args.mpi_ranks), *command]
    result, waves = {}, {}
    before = fd.geometry(inputs["restart"].read_bytes())
    for name, force in (("forces", "T"), ("electronic", "F")):
        work = args.output / name
        work.mkdir()
        (work / "si2.cntl").write_text(template.replace("!PSIDYN", f"!PSIDYN FORCE={force}"))
        (work / "model.fun").symlink_to(inputs["model"])
        shutil.copy2(inputs["structure"], work / "si2.strc")
        shutil.copy2(inputs["restart"], work / "si2.rstrt")
        fd.write_json(work / "inputs.json", {key: digest(work / key)
                      for key in ("si2.cntl", "si2.strc", "si2.rstrt", "model.fun")})
        with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
            subprocess.run(command, cwd=work, env=env, stdout=out, stderr=err,
                           timeout=1800, check=True)
        text = (work / "si2.prot").read_text()
        trace = stationarity.records(text)
        stationarity.validate(trace)
        if len(trace) != settings.block_steps:
            raise ValueError("Unexpected number of electronic steps")
        energies = fd.force_records(text, before["natom"]) if force == "T" else energy_records(text)
        if [row["step"] for row in energies] != [row["step"] for row in trace]:
            raise ValueError("Mismatched energy and electronic steps")
        result[name] = {"trace": trace, "energies": energies, "bands": stationarity.band_records(text),
                        "scalars": scalars(text, len(trace))}
        restart = (work / "si2.rstrt").read_bytes()
        fd.check_geometry(before, fd.geometry(restart))
        waves[name] = wave_records(restart)
        trajectory = work / "si2_f.tra"
        if force == "T" and (not trajectory.exists() or trajectory.stat().st_size == 0):
            raise ValueError("Missing requested force trajectory")
        if force == "F" and trajectory.exists() and trajectory.stat().st_size > 0:
            raise ValueError("Incomplete forces leaked into the trajectory")
        fd.write_json(work / "results.json", result[name])
        print(f"{name}: completed", flush=True)
    differences = compare(result["forces"], result["electronic"], args.tolerance)
    identical = waves["forces"] == waves["electronic"]
    if not identical:
        raise ValueError("Electronic restart records differ between force modes")
    for key, path in inputs.items():
        if digest(path) != provenance[key]["sha256"]:
            raise ValueError(f"Input changed during verification: {key}")
    summary = {"scope": "force-omission parity, not stationarity or force accuracy",
               "max_abs_difference": differences, "electronic_restart_bitwise_identical": identical,
               "tolerance": args.tolerance}
    fd.write_json(args.output / "results.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
