import copy
import unittest
from unittest.mock import patch

import parallel_grid as grid


class ParallelGridTest(unittest.TestCase):
    def test_rank_accounting_requires_every_rank_and_step(self):
        label = "GRID POINTS PER MPI RANK"
        self.assertEqual(grid.counters(label + " 3 2 0\n", label, 1, 3), [[3, 2, 0]])
        for text in ("", label + " 3 2\n", label + " 3 2 -1\n",
                     (label + " 3 2 0\n")*2):
            with self.assertRaises(ValueError):
                grid.counters(text, label, 1, 3)

    @patch.object(grid.modes, "wave_difference", return_value=0.)
    def test_comparison_rejects_missing_or_changed_diagnostics(self, unused):
        data = dict(scalars={"MODEL XC ENERGY": [-1.]},
                    energies=[dict(step=1, energy=-2.)], waves=[],
                    forces=[dict(forces=[(0., 0., 0.)])], stress=[(0., 0., 0.)]*3)
        self.assertEqual(max(grid.compare(data, data, "stress").values()), 0.)
        changes = [lambda d: d["scalars"].clear(),
                   lambda d: d["scalars"]["MODEL XC ENERGY"].clear(),
                   lambda d: d["energies"][0].update(energy=0.),
                   lambda d: d["forces"].clear(),
                   lambda d: d["forces"][0].update(forces=[(1., 0., 0.)]),
                   lambda d: d["stress"].clear(),
                   lambda d: d["stress"].__setitem__(0, (1., 0., 0.))]
        for change in changes:
            other = copy.deepcopy(data)
            change(other)
            with self.assertRaises(ValueError):
                grid.compare(data, other, "stress")


if __name__ == "__main__":
    unittest.main()
