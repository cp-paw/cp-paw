import unittest

import stationarity


def protocol(rms="2.1D-8", maximum="6.5E-8", commutator="1.5E-7", step=300):
    return f"""SKALA ELECTRONIC STATIONARITY STEP {step}
SKALA OCCUPIED RESIDUAL RMS {rms}
SKALA OCCUPIED RESIDUAL MAX {maximum}
SKALA OCCUPATION COMMUTATOR MAX {commutator}
SKALA SCF OVERLAP ERROR 2.5E-15
SKALA HAMILTONIAN HERMITICITY 4.2E-16
PROGRAM FINISHED
"""


class StationarityTest(unittest.TestCase):
    def test_warm_reference_meets_explicit_limits(self):
        result = stationarity.validate(stationarity.records(protocol()),
                                       residual=1e-6, commutator=1e-6)
        self.assertEqual(result["scope"], "stationarity")
        self.assertFalse(result["ground_state_minimum_certified"])

    def test_cold_snapshot_is_not_certified(self):
        data = stationarity.records(protocol("1.7", "1.9", "0.31"))
        self.assertEqual(stationarity.validate(data)["scope"], "diagnostic consistency")
        with self.assertRaisesRegex(ValueError, "maximum"):
            stationarity.validate(data, residual=1e-6, commutator=1e-6)

    def test_stationary_subspace_is_not_enough(self):
        with self.assertRaisesRegex(ValueError, "commutator"):
            stationarity.validate(stationarity.records(protocol(commutator="0.1")),
                                  residual=1e-6, commutator=1e-6)

    def test_invalid_values_and_incomplete_diagnostics(self):
        for raw in ("NaN", "Inf", "1e999", "-0.1", "1.2.3", "*****"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                stationarity.records(protocol(rms=raw))
        for text in (protocol().replace("PROGRAM FINISHED", ""),
                     protocol().replace("SKALA OCCUPIED RESIDUAL RMS 2.1D-8\n", ""),
                     "PROGRAM FINISHED"):
            with self.assertRaises(ValueError):
                stationarity.records(text)

    def test_duplicate_or_unscoped_value_fails(self):
        duplicate = "SKALA OCCUPIED RESIDUAL RMS 2.1D-8\n"
        for text in (duplicate + protocol(), protocol() + duplicate):
            with self.assertRaises(ValueError):
                stationarity.records(text)

    def test_every_step_checks_metric_and_finiteness(self):
        bad = protocol().replace("2.5E-15", "1E-2")
        data = stationarity.records(bad + protocol(step=301))
        with self.assertRaisesRegex(ValueError, "overlap"):
            stationarity.validate(data)
        with self.assertRaises(ValueError):
            stationarity.records(protocol(rms="NaN") + protocol(step=301))

    def test_final_window_is_explicit(self):
        data = stationarity.records(protocol("0.1", "0.2") + protocol(step=301))
        stationarity.validate(data, residual=1e-6, commutator=1e-6, last=1)
        with self.assertRaises(ValueError):
            stationarity.validate(data, residual=1e-6, commutator=1e-6, last=2)
        for count in (0, 3):
            with self.assertRaises(ValueError):
                stationarity.validate(data, last=count)

    def test_limits_must_be_complete_positive_and_finite(self):
        data = stationarity.records(protocol())
        for options in ({"residual": 1e-6}, {"commutator": 1e-6},
                        {"overlap": float("nan")}, {"hermiticity": 0.0},
                        {"residual": float("inf"), "commutator": 1e-6}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                stationarity.validate(data, **options)

    def test_rms_cannot_exceed_maximum(self):
        with self.assertRaisesRegex(ValueError, "RMS"):
            stationarity.records(protocol(rms="1.0"))


if __name__ == "__main__":
    unittest.main()
