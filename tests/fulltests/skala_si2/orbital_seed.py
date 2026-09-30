#!/usr/bin/env python3
"""Create a nonstationary periodic Si2 restart for orbital-derivative tests."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess

from cache_parity import digest


def structure(kmesh):
    if len(kmesh) != 3 or any(type(x) is not int or x < 1 for x in kmesh):
        raise ValueError("K mesh requires three positive integers")
    text = Path(__file__).with_name("skala_si2.strc").read_text()
    if text.count("DIV=2 2 2") != 1:
        raise ValueError("Si2 structure template changed")
    return text.replace("DIV=2 2 2", "DIV=" + " ".join(map(str, kmesh)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--kmesh", type=int, nargs=3, default=(3, 1, 1))
    parser.add_argument("--timeout", type=float, default=600)
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("Timeout must be finite and positive")
    try:
        text = structure(args.kmesh)
    except ValueError as error:
        parser.error(str(error))
    binary = args.executable.resolve(strict=True)
    work = args.output.resolve()
    work.mkdir(parents=True, exist_ok=False)
    (work / "si2.strc").write_text(text)
    (work / "si2.cntl").write_text("""!CONTROL
 !GENERIC TRACE=F DT=0.000001 NSTEP=1 NWRITE=1 START=T
          RSTRTTYPE='STATIC' AUTOCONV=1000 !END
 !DFT TYPE=10 !END
 !FOURIER EPWPSI=20 CDUAL=2 !END
 !CELL MOVE=F FRIC=0.0 M=1.E30 !END
 !PSIDYN MPSI=1000 FRIC=0.005 !END
!END
!EOB
""")
    env = dict(os.environ, CPPAW_GPU_MODE="off")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.pop("CPPAW_SKALA_ORBITAL_ROTATION", None)
    paths = (binary, work / "si2.cntl", work / "si2.strc")
    provenance = {"scope": "cold orbital probe, not a converged ground state",
                  "files": {str(path): digest(path) for path in paths},
                  "kmesh": args.kmesh,
                  "environment": {key: value for key, value in env.items()
                                  if key.startswith("CPPAW_") or key in
                                  ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")}}
    record = work / "provenance.json"
    record.write_text(json.dumps(provenance, indent=2) + "\n")
    with (work / "stdout.log").open("w") as out, (work / "stderr.log").open("w") as err:
        subprocess.run([str(binary), "si2.cntl"], cwd=work, env=env,
                       stdout=out, stderr=err, check=True, timeout=args.timeout)
    if "PROGRAM FINISHED" not in (work / "si2.prot").read_text():
        raise ValueError("Seed calculation did not finish normally")
    provenance["files"][str(work / "si2.rstrt")] = digest(work / "si2.rstrt")
    record.write_text(json.dumps(provenance, indent=2) + "\n")
    print(work)


if __name__ == "__main__":
    main()
