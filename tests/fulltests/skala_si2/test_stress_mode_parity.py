import copy
import unittest

import stress_mode_parity as modes
from test_force_stationary_fd import protocol


def stress_protocol(step=1):
    return protocol(step).replace("PROGRAM FINISHED", """SKALA TOTAL STRESS DIAGNOSTIC
TOTAL D E / D STRAIN 1.0 0.1 0.2
TOTAL D E / D STRAIN 0.1 2.0 0.3
TOTAL D E / D STRAIN 0.2 0.3 3.0
PROGRAM FINISHED""")


class StressModeTest(unittest.TestCase):
    def test_parser_associates_every_tensor_with_its_step(self):
        rows = modes.stress_records(stress_protocol() + stress_protocol(2))
        self.assertEqual([row["step"] for row in rows], [1, 2])
        self.assertEqual(rows[0]["tensor"], [[1., .1, .2], [.1, 2., .3], [.2, .3, 3.]])
        self.assertEqual(modes.stress_difference(rows, rows), 0.)
        changed = copy.deepcopy(rows)
        changed[1]["tensor"][2][1] += 0.01
        self.assertAlmostEqual(modes.stress_difference(rows, changed), .01)

    def test_parser_rejects_missing_duplicate_unscoped_and_nonfinite(self):
        text = stress_protocol()
        row = "TOTAL D E / D STRAIN 0.2 0.3 3.0\n"
        bad = [protocol(), text.replace(row, ""), text.replace(row, row*2),
               text.replace("SKALA TOTAL STRESS DIAGNOSTIC\n", ""),
               text.replace("TOTAL D E / D STRAIN 1.0", "TOTAL D E / D STRAIN nan"),
               text + "SKALA TOTAL STRESS DIAGNOSTIC\n",
               text + protocol(2), protocol() + stress_protocol(2),
               "SKALA TOTAL STRESS DIAGNOSTIC\n" + text]
        for value in bad:
            with self.subTest(value=value), self.assertRaises(ValueError):
                modes.stress_records(value)

    def test_comparison_rejects_changed_steps_shapes_and_nonfinite(self):
        original = modes.stress_records(stress_protocol())
        for change in (lambda d: d[0].update(step=2),
                       lambda d: d[0]["tensor"].pop(),
                       lambda d: d[0]["tensor"][0].pop(),
                       lambda d: d[0]["tensor"][0].__setitem__(0, float("nan"))):
            changed = copy.deepcopy(original)
            change(changed)
            with self.assertRaises(ValueError):
                modes.stress_difference(original, changed)
        with self.assertRaises(ValueError):
            modes.stress_difference([], [])


if __name__ == "__main__":
    unittest.main()
