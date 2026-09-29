import re
from types import SimpleNamespace
import unittest

import run


class CrystalDriverTest(unittest.TestCase):
    def test_reference_structures_and_real_arrays(self):
        for name, count in [("CO2", 12), ("NH3", 16), ("urea", 16)]:
            with self.subTest(crystal=name):
                lattice, atoms = run.read_structure(run.HERE / "structures" / f"{name}-solid.xyz")
                self.assertEqual(len(atoms), count)
                text = run.structure_text(lattice, atoms, 2)
                self.assertEqual(text.count(" !ATOM "), count)
                for element, _ in atoms:
                    self.assertIn(f"!SPECIES NAME='{element}_'", text)
                # CP-PAW's input parser rejects mixed integer/real arrays.
                for values in re.findall(r"\b(?:T|R)=([^!]+)!END", text):
                    self.assertTrue(all("e" in value for value in values.split()))

    def test_diagnostics_reject_missing_or_nonfinite_values(self):
        self.assertEqual(run.last_value("LABEL 1.23D-12\n", "LABEL"), 1.23e-12)
        for value in ["", "LABEL NaN\n", "LABEL 1e999\n"]:
            with self.assertRaises(ValueError):
                run.last_value(value, "LABEL")

    def test_probe_time_step_has_real_mantissa(self):
        args = SimpleNamespace(skala_steps=1, pbe_steps=180, device="CUDA",
                               radial=96, angular=17, image_shells=2, orientations=3,
                               dt=5.0, cutoff=40.0)
        for skala in [False, True]:
            text = run.control(args, skala)
            dt = re.search(r"\bDT=(\S+)", text).group(1)
            self.assertIn(".", dt.split("e")[0])
            if skala:
                self.assertIn("IMAGESHELLS=2", text)
                self.assertIn("LEBEDEVORIENTATIONS=3", text)

    def test_independent_counts_use_orthogonalizer_tolerance(self):
        record = {"OCCUPATION ELECTRONS": 8., "TRACE MINUS OCCUPATIONS": -9.4e-9,
                  "PS GRID MINUS TRACE": 2.e-14, "COMPOSITE MINUS TRACE": 0.1}
        # Atom-grid convergence is recorded separately and cannot be repaired
        # by passing the independent overlap/native-grid checks.
        run.check_electron_counts(record)
        for label in ["TRACE MINUS OCCUPATIONS", "PS GRID MINUS TRACE"]:
            with self.subTest(label=label), self.assertRaises(ValueError):
                run.check_electron_counts(dict(record, **{label: 1.e-4}))

    def test_higher_angular_request_is_preserved(self):
        args = SimpleNamespace(skala_steps=1, pbe_steps=180, device="CUDA",
                               radial=200, angular=64, image_shells=1, orientations=1,
                               dt=5.0, cutoff=40.0)
        # The Fortran grid library, not the input writer, rounds 64 up to 65.
        self.assertIn("LEBEDEVEXACTNESS=64", run.control(args, True))


if __name__ == "__main__":
    unittest.main()
