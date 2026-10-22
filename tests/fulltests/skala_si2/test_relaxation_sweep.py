import copy
import unittest

from relaxation_sweep import rank_settings


class RelaxationSweepTest(unittest.TestCase):
    def test_only_stationary_comparable_endpoints_can_be_selected(self):
        base = dict(name="a", elapsed_s=2., validated_steps=25,
                    result=dict(stationary=True, occupations=[1, 1], final=dict(energy=-2.)))
        faster = copy.deepcopy(base)
        faster.update(name="b", elapsed_s=1.)
        self.assertEqual(rank_settings([base, faster], 1e-7)["recommended"], "b")
        faster["validated_steps"] = 50
        self.assertEqual(rank_settings([base, faster], 1e-7)["recommended"], "a")
        faster["validated_steps"] = 25
        for update in (dict(stationary=False), dict(failure="failed residual"),
                       dict(final=dict(energy=-1.)), dict(final=dict(energy=float("nan"))),
                       dict(occupations=[1, 0])):
            other = copy.deepcopy(faster)
            other["result"].update(update)
            choice = rank_settings([base, other], 1e-7)["recommended"]
            self.assertIn(choice, ("a", None))
        self.assertIsNone(rank_settings([], 1e-7)["recommended"])


if __name__ == "__main__":
    unittest.main()
