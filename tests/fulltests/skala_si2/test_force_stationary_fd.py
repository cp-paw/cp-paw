import copy
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from displace_restart import displace, read_records, write_records
import force_stationary_fd as force_fd


def restart():
    def header(name):
        return struct.pack("<i", 4) + name.encode().ljust(32, b" ")
    xyz = struct.pack("<6d", 0, 0, 0, 2, 3, 4)
    cell = struct.pack("<27d", *([1, 0, 0, 0, 1, 0, 0, 0, 1] * 3))
    return write_records([header("CELL"), cell, header("ATOMS"), struct.pack("<i", 2),
                          b"names", xyz, xyz, b"unchanged electronic data"])


def protocol(step=1, residual="1e-8", energy="-2.0"):
    return f"""SKALA ELECTRONIC STATIONARITY STEP {step}
SKALA OCCUPIED RESIDUAL RMS {residual}
SKALA OCCUPIED RESIDUAL MAX {residual}
SKALA OCCUPATION COMMUTATOR MAX 1e-9
SKALA SCF OVERLAP ERROR 1e-12
SKALA HAMILTONIAN HERMITICITY 1e-14
SKALA BAND RESIDUAL 1 1 1 1.0 {residual} -1.0 1e-9
SKALA TOTAL FORCE DIAGNOSTIC
============================
TOTAL ENERGY {energy}
ATOM 1 0.1 0.2 0.3
ATOM 2 -0.1 -0.2 -0.3
NET FORCE 0 0 0
PROGRAM FINISHED
"""


class StationaryForceTest(unittest.TestCase):
    def test_relaxation_retains_failed_blocks_and_requires_a_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = {key: root / key for key in ("restart", "model", "structure", "executable")}
            for path in inputs.values():
                path.write_bytes(b"fixture")
            inputs["restart"].write_bytes(restart())
            args = SimpleNamespace(output=root, atom=2, axis=1, dt=5, block_steps=1,
                                   max_blocks=2, last=1, residual_tolerance=1e-6,
                                   commutator_tolerance=1e-6, device="CPU", radial_points=96,
                                   lebedev_exactness=17, cutoff=20, mass=25,
                                   mass_g2=0.3166286988823056, friction=0.05, timeout=1)
            def run(command, cwd, **kwargs):
                residual = "0.1" if cwd.name == "block-1" else "1e-8"
                (cwd / "si2.prot").write_text(protocol(residual=residual))
            with patch.object(force_fd.subprocess, "run", side_effect=run):
                _, result = force_fd.run_leg("center", 0, args, inputs, {})
            self.assertTrue(result["stationary"])
            self.assertEqual([block["stationary"] for block in result["blocks"]], [False, True])
            failed = json.loads((root / "center/block-1/results.json").read_text())
            self.assertFalse(failed["stationary"])
            self.assertIn("stationarity_failure", failed)
            args.max_blocks = 1
            with patch.object(force_fd.subprocess, "run", side_effect=run):
                _, result = force_fd.run_leg("unconverged", 0, args, inputs, {})
            self.assertFalse(result["stationary"])
            self.assertIn("failure", result)
            args.rigid_translation = True
            with patch.object(force_fd.subprocess, "run", side_effect=run):
                force_fd.run_leg("translated", 0.003, args, inputs, {})
            translated = force_fd.geometry((root / "translated/block-1/si2.rstrt").read_bytes())
            self.assertEqual(translated["positions"], [[0.003, 0, 0, 2.003, 3, 4]] * 2)

    def test_control_has_typed_real_parameters_and_fixed_geometry(self):
        args = SimpleNamespace(dt=5, block_steps=40, device="CPU", radial_points=96,
                               lebedev_exactness=17, cutoff=20, mass=25,
                               mass_g2=0.3166286988823056, friction=0.05)
        text = force_fd.control(args)
        self.assertIn("DT=5.0000000000000000e+00", text)
        self.assertIn("!CELL MOVE=F FRIC=0.0", text)
        self.assertNotIn("!RDYN", text)
        self.assertNotIn("!MERMIN", text)
        self.assertIn("SAFEORTHO=T", text)
        self.assertNotIn("ORTHOTOL", text)
        args.orthogonality_tolerance = 1e-12
        self.assertIn("ORTHOTOL=9.9999999999999998e-13", force_fd.control(args))

    def test_invalid_orthogonality_tolerances_fail_before_creating_output(self):
        arguments = ["force_stationary_fd.py"]
        for key in ("executable", "model", "restart", "structure", "output"):
            arguments.extend(["--" + key, "unused"])
        for value in ("0", "1e-15", "1e-7", "nan", "inf"):
            with self.subTest(value=value), redirect_stderr(io.StringIO()), \
                    patch.object(sys, "argv", arguments + ["--orthogonality-tolerance", value]), \
                    self.assertRaises(SystemExit) as error:
                force_fd.main()
            self.assertEqual(error.exception.code, 2)

    def test_displacement_preserves_other_records_and_both_time_levels(self):
        initial = restart()
        changed = displace(initial, 2, 1, 0.003)
        before, after = read_records(initial), read_records(changed)
        for index, (old, new) in enumerate(zip(before, after)):
            if index in (5, 6):
                self.assertEqual(struct.unpack("<6d", new), (0, 0, 0, 2.003, 3, 4))
            else:
                self.assertEqual(old, new)
        self.assertEqual(force_fd.geometry(initial)["natom"], 2)
        with self.assertRaisesRegex(ValueError, "moved"):
            force_fd.check_geometry(force_fd.geometry(initial), force_fd.geometry(changed))

    def test_restart_rejects_bad_geometry_and_displacements(self):
        for atom, axis, delta in ((0, 1, 0), (3, 1, 0), (1, 4, 0), (1, 1, float("nan"))):
            with self.assertRaises(ValueError):
                displace(restart(), atom, axis, delta)
        for data in (b"", restart()[:-1]):
            with self.assertRaises(ValueError):
                force_fd.geometry(data)

    def test_rigid_translation_preserves_cell_and_electronic_records(self):
        initial = restart()
        for axis in (1, 2, 3):
            with self.subTest(axis=axis):
                changed = force_fd.displace_geometry(initial, 2, axis, -0.003, True)
                before, after = read_records(initial), read_records(changed)
                for index, (old, new) in enumerate(zip(before, after)):
                    if index in (5, 6):
                        expected = list(struct.unpack("<6d", old))
                        for atom in (0, 1):
                            expected[3 * atom + axis - 1] -= 0.003
                        self.assertEqual(struct.unpack("<6d", new), tuple(expected))
                    else:
                        self.assertEqual(old, new)
        self.assertEqual(force_fd.displace_geometry(initial, 2, 1, 0.003),
                         displace(initial, 2, 1, 0.003))
        for axis, delta in ((0, 0), (4, 0), (1, float("nan")), (1, float("inf"))):
            with self.subTest(axis=axis, delta=delta), self.assertRaises(ValueError):
                force_fd.displace_geometry(initial, 2, axis, delta, True)

    def test_rigid_translation_compares_the_sum_of_all_forces(self):
        center = force_fd.analyze(protocol(), 2, 1, 1, 1e-6, 1e-6)
        center["final"]["forces"] = [[0.1, 0.2, 0.3], [-0.06, -0.1, -0.2]]
        minus, plus = copy.deepcopy(center), copy.deepcopy(center)
        minus["final"]["energy"] = -1.99996
        plus["final"]["energy"] = -2.00004
        row = force_fd.compare(center, minus, plus, 0.001, 2, 1, 1e-10,
                               rigid_translation=True)
        self.assertAlmostEqual(row["analytic_hartree_per_bohr"], 0.04)
        self.assertAlmostEqual(row["finite_difference_hartree_per_bohr"], 0.04)
        self.assertTrue(row["passed"])
        self.assertFalse(force_fd.compare(center, minus, plus, 0.001, 2, 1, 1e-10)["passed"])
        for atom, axis in ((0, 1), (3, 1), (1, 0), (1, 4)):
            with self.assertRaises(ValueError):
                force_fd.force_component(center, atom, axis)
        plus["stationary"] = False
        with self.assertRaisesRegex(ValueError, "stationarity"):
            force_fd.compare(center, minus, plus, 0.001, 2, 1, rigid_translation=True)

    def test_parser_uses_last_force_and_energy_not_first(self):
        text = protocol() + protocol(step=2, energy="-2.1D0")
        result = force_fd.analyze(text, 2, 2, 2, 1e-6, 1e-6)
        self.assertTrue(result["stationary"])
        self.assertEqual(result["final"]["energy"], -2.1)
        self.assertEqual(result["final"]["step"], 2)
        self.assertAlmostEqual(result["final_energy_span_hartree"], 0.1)

    def test_force_energy_retains_binary64_roundtrip_digits(self):
        value = -46.938882201616096
        text = protocol(energy=f"{value:.16E}")
        result = force_fd.force_records(text, 2)
        self.assertEqual(result[0]["energy"], value)
        self.assertNotEqual(float(f"{value:.14E}"), value)

    def test_nonstationarity_is_retained_and_blocks_comparison(self):
        text = protocol(residual="0.1") + protocol(step=2)
        failed = force_fd.analyze(text, 2, 2, 2, 1e-6, 1e-6)
        self.assertFalse(failed["stationary"])
        self.assertIn("stationarity_failure", failed)
        with self.assertRaisesRegex(ValueError, "stationarity"):
            force_fd.compare(failed, failed, failed, 0.001, 2, 1)
        self.assertTrue(force_fd.analyze(text, 2, 2, 1, 1e-6, 1e-6)["stationary"])

    def test_missing_duplicate_nonfinite_and_unscoped_reports_fail(self):
        text = protocol()
        variants = [text.replace("TOTAL ENERGY -2.0", "TOTAL ENERGY NaN"),
                    text.replace("ATOM 2 -0.1", "ATOM 2 NaN"),
                    text.replace("ATOM 2 -0.1", "ATOM 1 -0.1"),
                    text.replace("TOTAL ENERGY -2.0", "TOTAL ENERGY -2.0\nTOTAL ENERGY -2.0"),
                    text.replace("NET FORCE 0 0 0", "NET FORCE 1 0 0"),
                    text.replace("NET FORCE 0 0 0", ""),
                    text.replace("SKALA ELECTRONIC STATIONARITY STEP 1", ""),
                    text.replace("SKALA TOTAL FORCE DIAGNOSTIC", ""),
                    text.replace("ATOM 2 -0.1 -0.2 -0.3", ""),
                    text + text[text.index("SKALA TOTAL FORCE DIAGNOSTIC"):]]
        for variant in variants:
            with self.subTest(text=variant), self.assertRaises(ValueError):
                force_fd.force_records(variant, 2)

    def test_step_band_and_completion_checks_are_mandatory(self):
        for text in (protocol().replace("PROGRAM FINISHED", ""),
                     protocol() + protocol(step=1), protocol() + protocol(step=3),
                     protocol() + protocol(step=2).replace("1 1 1 1.0", "1 1 1 0.5"),
                     protocol() + protocol(step=2).replace("1e-12", "1e-3")):
            with self.assertRaises(ValueError):
                force_fd.analyze(text, 2, 2, 1, 1e-6, 1e-6)

    def test_force_sign_and_explicit_acceptance(self):
        center = force_fd.analyze(protocol(), 2, 1, 1, 1e-6, 1e-6)
        minus, plus = copy.deepcopy(center), copy.deepcopy(center)
        minus["final"]["energy"] = -2.0001
        plus["final"]["energy"] = -1.9999
        row = force_fd.compare(center, minus, plus, 0.001, 2, 1)
        self.assertAlmostEqual(row["finite_difference_hartree_per_bohr"], -0.1)
        self.assertIsNone(row["passed"])
        self.assertTrue(force_fd.compare(center, minus, plus, 0.001, 2, 1, 1e-10)["passed"])
        plus["occupations"] = []
        with self.assertRaisesRegex(ValueError, "occupations"):
            force_fd.compare(center, minus, plus, 0.001, 2, 1)


if __name__ == "__main__":
    unittest.main()
