from types import SimpleNamespace
import unittest

import run


class NativeScanDriverTest(unittest.TestCase):
    def test_structure_controls_setups_without_an_ignored_override(self):
        args = SimpleNamespace(electronic_only=False, stress=True, dt=1e-6, steps=1,
            radial_points=200, lebedev_exactness=53, cutoff=20., density_dual=4.,
            mass=25., mass_g2=0.3166286988823056, friction=0.4)
        text = run.control(args)
        self.assertNotIn("!FILES", text)
        self.assertNotIn("stp.cntl", text)
        self.assertNotIn("PARMS_STP", text)
        self.assertIn("!SCAN RADIALPOINTS=200 LEBEDEVEXACTNESS=53", text)


if __name__ == "__main__":
    unittest.main()
