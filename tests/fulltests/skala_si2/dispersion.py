#!/usr/bin/env python3
"""Periodic Si2 D3(BJ) increments against Grimme's original dftd3 executable."""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np


STRUCTURE = """!STRUCTURE
 !GENERIC LUNIT=1.0 STPVERSION='2.0.0' !END
 !KPOINTS DIV=1 1 1 !END
 !OCCUPATIONS EMPTY=5 NSPIN=1 !END
 !LATTICE T=0.0 5.13 5.13  5.13 0.0 5.13  5.13 5.13 0.0 !END
 !SPECIES NAME='Si' ZV=4. M=5. NPRO=2 2 1 LRHOX=2 ID='SI_.75_6.0' !END
 !ATOM NAME='Si_1' R=0.0 0.0 0.0 !END
 !ATOM NAME='Si_2' R=2.68 2.38 2.51 !END
!END
!EOB
"""

# Original Bonn dftd3 (June 2016, not s-dftd3):
# dftd3 POSCAR -pbc -bj (or -zero) -func b3-lyp (or pbe) -cutoff 75.58905 -cnthr 40 -grad
# POSCAR used Cartesian coordinates as above and its native 0.52917726 scale.
# Energies are printed to 8 decimals. Rows of dftd3_cellgradient correspond to
# lattice vectors. Contract with the lattice to obtain dE/d(strain), in Ha.
E_BJ = -0.04348184
E_PBE_BJ = -0.02313829
F1 = np.array([-0.00037578519184579, 0.00058642392419702, 0.00021335394488323])
CELL = np.array([[0, 5.13, 5.13], [5.13, 0, 5.13], [5.13, 5.13, 0]])
STRAIN = CELL.T @ np.array([
    [-0.0041496480338034, 0.0041636905215205, 0.0041413154757329],
    [0.0043234951646394, -0.0042853553696135, 0.0044368274546618],
    [0.0042430399419265, 0.0043787472777366, -0.0042629803238258]])
ZERO_F1 = np.array([-0.000033670094314918, 0.000053791194651845, 0.000017751937332810])
ZERO_STRAIN = CELL.T @ np.array([
    [-0.00084176989685273, 0.00085496548618414, 0.00082992608846099],
    [0.0010176748144455, -0.00097295972291768, 0.0011238252427873],
    [0.00092073087449150, 0.0010519207005564, -0.00094792032519453]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    args = parser.parse_args()
    cmd = [str(args.binary.resolve(strict=True)), "si2.cntl"]
    if args.mpi_ranks > 1:
        cmd = ["mpirun", "-np", str(args.mpi_ranks)] + cmd
    model = args.model.resolve(strict=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="paw-d3bj-") as tmp:
        root = Path(tmp)

        def run(name, vdw, apply=None, seed=False, scan=False, structure=STRUCTURE):
            work = root / name
            work.mkdir()
            skala = "" if apply is None else (
                f"!SKALA MODEL='model.fun' DEVICE='{args.device.upper()}' RADIALPOINTS=32 LEBEDEVEXACTNESS=17 "
                f"APPLY={apply} CHECK=T !END")
            if scan:
                skala = "!SCAN RADIALPOINTS=32 LEBEDEVEXACTNESS=17 CHECK=T !END"
            (work / "si2.strc").write_text(structure)
            (work / "model.fun").symlink_to(model)
            if not seed:
                shutil.copyfile(root / "seed/si2.rstrt", work / "si2.rstrt")
            (work / "si2.cntl").write_text(f"""!CONTROL
 !GENERIC TRACE=F DT=0.000001 NSTEP=1 NWRITE=1 START={'T' if seed else 'F'} RSTRTTYPE='STATIC' !END
 !DFT TYPE=10 VDW={vdw} {skala} !END
 !FOURIER EPWPSI=20 CDUAL=2 !END
 !CELL MOVE=F FRIC=0. M=1.E30 !END
 !PSIDYN FORCE=T STRESS=T MPSI=1000 FRIC=0.005 !END
!END
!EOB
""")
            p = subprocess.run(cmd, cwd=work, env=env, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=360)
            prot = work / "si2.prot"
            text = prot.read_text() if prot.exists() else ""
            if p.returncode or "PROGRAM FINISHED" not in text:
                raise RuntimeError(f"{name}: {p.stdout[-3000:]}\n{text[-2000:]}")
            return text

        run("seed", "F", seed=True)
        pbe = run("pbe", "BJ")
        def term(text, method, damping="BJ"):
            m = re.search(rf"^D3\({damping}\) {method} ATM=F PERIODIC=T ENERGY \[HA\]\s+(\S+)", text, re.M | re.I)
            if not m:
                raise ValueError(f"Missing {method} D3({damping}) term")
            return float(m[1])
        np.testing.assert_allclose(term(pbe, "pbe"), E_PBE_BJ, rtol=0, atol=1e-8)
        # Legacy PBE before the BJ integration (protocol printed to 7 decimals).
        # The legacy implementation also fits atomic charges. A single Gaussian
        # per site avoids a singular fit on this deliberately coarse FFT grid.
        legacy = run("pbe-legacy", "T", structure=STRUCTURE.replace(
            " !KPOINTS", " !ISOLATE DECOUPLE=F NF=1 !END\n !KPOINTS"))
        energy = re.search(r"^VAN DER WAALS ENERGY\s*:\s*([-+.\dEe]+)", legacy, re.M)
        assert energy and "D3(BJ)" not in legacy
        np.testing.assert_allclose(float(energy[1]), -0.0062751, rtol=0, atol=5e-8)
        # The 2016 dftd3 release lacks SCAN; use the standalone s-dftd3 1.7.0 result.
        scan = run("scan", "'BJ'", scan=True)
        np.testing.assert_allclose(term(scan, "scan"), -0.00728046904255, rtol=0, atol=1e-9)
        scan_zero = run("scan-zero", "T", scan=True)
        np.testing.assert_allclose(term(scan_zero, "scan", "0"), -0.00393670949126, rtol=0, atol=1e-9)
        for apply in ("T", "F"):
            off = run("off-" + apply, "F", apply)
            for mode in ("T", "BJ"):
                on = run("on-" + apply + "-" + mode, mode, apply)
                damping, er, fr, sr = (("BJ", E_BJ, F1, STRAIN) if mode == "BJ" else
                                       ("0", -0.01853348, ZERO_F1, ZERO_STRAIN))
                ed3 = term(on, "b3lyp", damping)
                def values(text):
                    e = float(re.search(r"^TOTAL ENERGY\s+(\S+)\s*$", text, re.M)[1])
                    def matrix(pattern):
                        return np.array([list(map(float, s.split())) for s in re.findall(pattern, text, re.M)])
                    return (e, matrix(r"^ATOM\s+\d+\s+([-+.\dEe ]+)$"),
                            matrix(r"^TOTAL D E / D STRAIN\s+([-+.\dEe ]+)$"))
                e0, f0, s0 = values(off)
                e1, f1, s1 = values(on)
                de, df, ds = e1-e0, f1-f0, s1-s0
                assert df.shape == (2,3) and ds.shape == (3,3)
                errors = abs(de-er), np.max(abs(df-[fr,-fr])), np.max(abs(ds-sr))
                print(f"APPLY={apply} VDW={mode}: delta E={de:.12g}, max reference errors E/F/dEdeps:",
                      *errors, flush=True)
                np.testing.assert_allclose(de, ed3, rtol=0, atol=2e-9)
                np.testing.assert_allclose(de, er, rtol=0, atol=2e-8 if mode == "T" else 1e-8)
                np.testing.assert_allclose(df, [fr,-fr], rtol=0, atol=1e-8)
                np.testing.assert_allclose(ds, sr, rtol=0, atol=5e-8)


if __name__ == "__main__":
    main()
