import copy
import unittest

import orbital_fd
import orbital_seed
from test_stationarity import protocol


def sample():
    return protocol() + """SKALA ORBITAL ROTATION INDICES 1 1 4 5
SKALA ORBITAL ROTATION MODE REAL
SKALA ORBITAL ROTATION ANGLE 0.02
SKALA ORBITAL ROTATION OCCUPATIONS 0.25 0
SKALA ORBITAL ROTATION HAMILTONIAN 0.4 0.2
SKALA ORBITAL ROTATION DERIVATIVE 0.2
SKALA ORBITAL ROTATION PACKED 1
SKALA TOTAL FORCE DIAGNOSTIC
TOTAL ENERGY -100.000000000000
NET FORCE 0 0 0
TOTAL ENERGY : -100.0000
"""


class OrbitalDerivativeTest(unittest.TestCase):
    def test_complex_seed_mesh(self):
        self.assertIn("DIV=3 1 1", orbital_seed.structure((3, 1, 1)))
        for mesh in ((0, 1, 1), (-1, 1, 1), (1, 1), (3.0, 1, 1)):
            with self.subTest(mesh=mesh), self.assertRaises(ValueError):
                orbital_seed.structure(mesh)

    def test_weighted_occupations_and_high_precision_energy(self):
        result = orbital_fd.diagnostics(sample())
        self.assertEqual(result["DERIVATIVE"], [0.2])
        self.assertEqual(result["ENERGY"], [-100])

    def test_imaginary_sign(self):
        text = sample().replace("MODE REAL", "MODE IMAG").replace("PACKED 1", "PACKED 0")
        text = text.replace("DERIVATIVE 0.2", "DERIVATIVE -0.1")
        self.assertEqual(orbital_fd.diagnostics(text)["DERIVATIVE"], [-0.1])
        with self.assertRaises(ValueError):
            orbital_fd.diagnostics(text.replace("PACKED 0", "PACKED 1"))

    def test_missing_duplicate_malformed_and_nonfinite(self):
        line = "SKALA ORBITAL ROTATION ANGLE 0.02\n"
        for text in (sample().replace(line, ""), sample() + line,
                     sample().replace("ANGLE 0.02", "ANGLE NaN"),
                     sample().replace("ANGLE 0.02", "ANGLE 1e999"),
                     sample().replace("INDICES 1 1 4 5", "INDICES 1 1 4.5 5"),
                     sample().replace("DERIVATIVE 0.2", "DERIVATIVE 0.1"),
                     sample().replace("OCCUPATIONS 0.25 0", "OCCUPATIONS 0 0"),
                     sample().replace("PROGRAM FINISHED", "")):
            with self.subTest(text=text), self.assertRaises(ValueError):
                orbital_fd.diagnostics(text)

    def test_centered_difference_with_explicit_limit(self):
        center = orbital_fd.diagnostics(sample())
        minus, plus = copy.deepcopy(center), copy.deepcopy(center)
        minus["ANGLE"], plus["ANGLE"] = [0.01], [0.03]
        minus["ENERGY"], plus["ENERGY"] = [-100.002], [-99.998]
        row = orbital_fd.compare(center, minus, plus, 0.01, 1e-10)
        self.assertTrue(row["passed"])
        self.assertIsNone(orbital_fd.compare(center, minus, plus, 0.01)["passed"])
        plus["ENERGY"] = [-99.997]
        self.assertFalse(orbital_fd.compare(center, minus, plus, 0.01, 1e-10)["passed"])

    def test_changed_occupations_and_angles_rejected(self):
        center = orbital_fd.diagnostics(sample())
        minus, plus = copy.deepcopy(center), copy.deepcopy(center)
        minus["ANGLE"], plus["ANGLE"] = [0.01], [0.03]
        for key, value in (("OCCUPATIONS", [0.2, 0]), ("INDICES", [2, 1, 4, 5]),
                           ("ANGLE", [0.02]), ("MODE", "IMAG"), ("PACKED", [0])):
            bad = copy.deepcopy(plus)
            bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                orbital_fd.compare(center, minus, bad, 0.01)
        for step in (0, -0.01, float("nan")):
            with self.assertRaises(ValueError):
                orbital_fd.compare(center, minus, plus, step)


if __name__ == "__main__":
    unittest.main()
