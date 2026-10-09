#!/usr/bin/env python3
"""Require one incremental build to relink after recompiling a generated source.

Run only in an idle build directory. The test touches its generated Fortran
source, never the repository source, and recompiles the affected module.
"""

import argparse
import os
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--make", default="make")
    args = parser.parse_args()
    build = args.build.resolve(strict=True)
    source = build / "paw_waves2.f90"
    obj = build / "paw_waves2.o"
    binary = build / "paw.x"
    for path in (source, obj, binary):
        if not path.is_file():
            parser.error(f"Build CP-PAW first: missing {path}")
    previous = obj.stat().st_mtime_ns
    time.sleep(1.1)
    os.utime(source, None)
    subprocess.run([args.make, "-s", "-j2", "executable"], cwd=build, check=True)
    if obj.stat().st_mtime_ns <= previous:
        raise RuntimeError("The changed generated source was not recompiled")
    if binary.stat().st_mtime_ns < obj.stat().st_mtime_ns:
        raise RuntimeError("One incremental invocation left a stale executable")
    linked = binary.stat().st_mtime_ns
    subprocess.run([args.make, "-s", "-j2", "executable"], cwd=build, check=True)
    if binary.stat().st_mtime_ns != linked:
        raise RuntimeError("An unchanged build unnecessarily relinked the executable")
    print("Incremental source compilation, relink and unchanged rebuild: passed")


if __name__ == "__main__":
    main()
