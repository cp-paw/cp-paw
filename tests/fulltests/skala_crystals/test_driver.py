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
                               radial=96, angular=17, image_shells=2, dt=5.0, cutoff=40.0)
        for skala in [False, True]:
            text = run.control(args, skala)
            dt = re.search(r"\bDT=(\S+)", text).group(1)
            self.assertIn(".", dt.split("e")[0])
            if skala:
                self.assertIn("IMAGESHELLS=2", text)


if __name__ == "__main__":
    unittest.main()
