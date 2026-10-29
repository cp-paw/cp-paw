from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import relax
import run


class NativeScanDriverTest(unittest.TestCase):
    def test_relaxation_checks_both_functionals_without_relaxing_the_gate(self):
        for functional in ("scan", "skala"):
            for stationary in (False, True):
                with self.subTest(functional=functional, stationary=stationary), \
                        tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    argv = ["relax.py", "--functional", functional, "--electronic-only",
                            "--output", str(root / "result")]
                    keys = ("executable", "restart", "structure")
                    if functional == "skala":
                        keys += ("model",)
                    for key in keys:
                        path = root / key
                        path.write_bytes(b"fixture")
                        argv.extend(["--" + key, str(path)])

                    def leg(name, delta, args, inputs, env, *, energy_only):
                        self.assertEqual((name, delta, args.functional), ("relax", 0., functional))
                        self.assertEqual(set(inputs), set(keys))
                        self.assertTrue(energy_only)
                        self.assertEqual((args.residual_tolerance, args.commutator_tolerance),
                                         (1e-7, 1e-7))
                        control = relax.fd.control(args, electronic_only=True)
                        self.assertIn("START=F", control)
                        self.assertIn("!SCAN" if functional == "scan" else "!SKALA MODEL=", control)
                        return name, {"blocks": [], "stationary": stationary}

                    with patch.object(sys, "argv", argv), \
                            patch.object(relax.fd, "run_leg", side_effect=leg), \
                            redirect_stdout(io.StringIO()):
                        self.assertEqual(relax.main(), 0 if stationary else 1)
                    result = json.loads((root / "result/results.json").read_text())
                    self.assertEqual(result["functional"], functional)
                    self.assertEqual(result["stationary"], stationary)
                    self.assertFalse(result["forces_evaluated"])

    def test_kernel_runs_on_every_nonempty_rank(self):
        text = "GRID POINTS PER MPI RANK 3 3 2 0\nMODEL CALLS PER MPI RANK 1 1 1 0\n"
        self.assertEqual(run.parallel_work(text, 1, 4),
                         {"rows": [[3, 3, 2, 0]], "calls": [[1, 1, 1, 0]]})
        for broken in (text.replace("1 1 1 0", "1 1 0 0"),
                       text.replace("1 1 1 0", "1 1 1 1"),
                       text.replace("3 3 2 0", "3 3 -2 0")):
            with self.subTest(text=broken), self.assertRaises(ValueError):
                run.parallel_work(broken, 1, 4)
        with self.assertRaises(ValueError):
            run.parallel_work(text, 2, 4)
        with self.assertRaises(ValueError):
            run.parallel_work(text, 1, 3)

    def test_structure_controls_setups_without_an_ignored_override(self):
        args = SimpleNamespace(electronic_only=False, stress=True, dt=1e-6, steps=1,
            radial_points=200, lebedev_exactness=53, cutoff=20., density_dual=4.,
            mass=25., mass_g2=0.3166286988823056, friction=0.4)
        text = run.control(args)
        self.assertNotIn("!FILES", text)
        self.assertNotIn("stp.cntl", text)
        self.assertNotIn("PARMS_STP", text)
        self.assertIn("!SCAN RADIALPOINTS=200 LEBEDEVEXACTNESS=53", text)
        self.assertIn("LEBEDEVORIENTATIONS=1", text)
        args.lebedev_orientations = 4
        self.assertIn("LEBEDEVORIENTATIONS=4", run.control(args))


if __name__ == "__main__":
    unittest.main()
