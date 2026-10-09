#!/usr/bin/env python3
"""A failed generator or compiler must not leave an up-to-date partial target."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    make = os.environ.get("MAKE") or shutil.which("gmake") or shutil.which("make")
    root = Path(__file__).resolve().parents[2]
    template = root / "src/Buildtools/Makefile.in"
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        source = work / "probe.f90pp"
        target = work / "probe.f90"
        source.write_text("program probe\nend program\n")
        processor = work / "preprocess.sh"
        processor.write_text("#!/bin/sh\nprintf 'partial output\\n'\nexit 1\n")
        processor.chmod(0o755)
        command = [make, "-s", "-f", str(template), "-o", "preprocess.sh",
                   "F90PP=./preprocess.sh", "CPPFLAGS=", "probe.f90"]
        failed = subprocess.run(command, cwd=work, capture_output=True, text=True)
        assert failed.returncode != 0, failed.stdout
        assert "Deleting file 'probe.f90'" in failed.stderr, failed.stderr
        assert not target.exists(), "Failed preprocessing left a partial generated source"
        processor.write_text("#!/bin/sh\ncat\n")
        subprocess.run(command, cwd=work, check=True)
        assert target.read_bytes() == source.read_bytes()
        previous = target.stat().st_mtime_ns
        subprocess.run(command, cwd=work, check=True)
        assert target.stat().st_mtime_ns == previous

        # Exercise the actual generated recursive makefile too, without a compiler.
        (work / "probe.mk").write_text("probe.o:\n\tprintf 'partial object' > $@\n\tfalse\n")
        subprocess.run([make, "-s", "-f", str(template), "-o", "probe.mk",
                        "PAWLIST=probe", "LIBLIST=", "TOOLS=",
                        "FC=unused", "FCFLAGS=", "big.mk"], cwd=work, check=True)
        failed = subprocess.run([make, "-s", "-f", "big.mk", "probe.o"], cwd=work,
                                capture_output=True, text=True)
        assert failed.returncode != 0, failed.stdout
        assert not (work / "probe.o").exists(), "Failed compilation left a partial object"
    print("Failed-build cleanup, recovery and unchanged rebuild: passed")


if __name__ == "__main__":
    main()
