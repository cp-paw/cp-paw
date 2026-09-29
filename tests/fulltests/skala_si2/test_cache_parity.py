import unittest

import cache_parity
import stationarity


def protocol(cached=False):
    lines = []
    for step in range(2):
        lines.append(f"{stationarity.HEADER} {step}")
        for label, (width, rows) in cache_parity.LAYOUT.items():
            for _ in range(rows // 2):
                lines.append(label + " " + " ".join(["0.0"] * width))
        lines.extend(["ATOM-GRID ROWS 100",
                      f"PARTITION CACHE HITS {100 if cached and step else 0}",
                      f"PARTITION CACHE MISSES {0 if cached and step else 100}",
                      f"PARTITION CACHE BYTES ALL RANKS {10000 if cached else 0}"])
    return "\n".join([*lines, "PROGRAM FINISHED"])


class CacheParityTests(unittest.TestCase):
    def test_equal(self):
        result = cache_parity.compare(cache_parity.diagnostics(protocol()),
                                      cache_parity.diagnostics(protocol(True)), 1e-10)
        self.assertEqual(result["rows_per_step"], 100)
        self.assertEqual(max(result["max_abs_difference"].values()), 0)

    def test_real_force_section(self):
        text = protocol().replace("TOTAL FORCE ATOM 0.0 0.0 0.0 0.0",
                                   "SKALA TOTAL FORCE DIAGNOSTIC\n"
                                   "============================\n"
                                   "ATOM 0.0 0.0 0.0 0.0\nNET FORCE 0 0 0")
        self.assertEqual(cache_parity.diagnostics(text), cache_parity.diagnostics(protocol()))

    def test_missing_or_duplicate_row(self):
        for text in (protocol().replace("MODEL XC ENERGY 0.0\n", "", 1),
                     protocol() + "\nMODEL XC ENERGY 0.0"):
            with self.assertRaises(ValueError):
                cache_parity.diagnostics(text)

    def test_hidden_nonfinite_or_malformed(self):
        for value in ("NaN", "Inf", "1e999", "oops", "0.0 trailing"):
            with self.assertRaises(ValueError):
                cache_parity.diagnostics(protocol().replace("MODEL XC ENERGY 0.0",
                                                           "MODEL XC ENERGY " + value, 1))

    def test_cache_must_really_hit(self):
        with self.assertRaises(ValueError):
            cache_parity.compare(cache_parity.diagnostics(protocol()),
                                 cache_parity.diagnostics(protocol()), 1e-10)

    def test_invalid_counters(self):
        for value in ("-1", "0.5"):
            with self.assertRaises(ValueError):
                cache_parity.diagnostics(protocol().replace("PARTITION CACHE HITS 0",
                                                           "PARTITION CACHE HITS " + value, 1))

    def test_stress_can_invalidate_geometry(self):
        changed = protocol(True).replace("PARTITION CACHE HITS 100", "PARTITION CACHE HITS 0")
        changed = changed.replace("PARTITION CACHE MISSES 0", "PARTITION CACHE MISSES 100")
        result = cache_parity.compare(cache_parity.diagnostics(protocol(), stress=True),
                                     cache_parity.diagnostics(changed, stress=True), 1e-10,
                                     require_reuse=False)
        self.assertIn("TOTAL D E / D STRAIN", result["max_abs_difference"])
        self.assertEqual(result["cache_hits_per_step"], [0, 0])

    def test_force_difference(self):
        changed = protocol(True).replace("SKALA FORCE ATOM 0.0 0.0 0.0 0.0",
                                          "SKALA FORCE ATOM 0.0 0.0 0.1 0.0", 1)
        with self.assertRaises(ValueError):
            cache_parity.compare(cache_parity.diagnostics(protocol()),
                                 cache_parity.diagnostics(changed), 1e-10)

    def test_total_force_has_separate_bound(self):
        off = cache_parity.diagnostics(protocol())
        cached = cache_parity.diagnostics(protocol(True))
        cached["TOTAL FORCE ATOM"][2] = (0., 0., 5e-9, 0.)
        cache_parity.compare(off, cached, 1e-10)
        with self.assertRaises(ValueError):
            cache_parity.compare(off, cached, 1e-10, total_force_tolerance=1e-10)
        cached["TOTAL FORCE ATOM"][2] = (0., 0., 2e-8, 0.)
        with self.assertRaises(ValueError):
            cache_parity.compare(off, cached, 1e-10)

    def test_unfinished(self):
        with self.assertRaises(ValueError):
            cache_parity.diagnostics(protocol().replace("PROGRAM FINISHED", ""))


if __name__ == "__main__":
    unittest.main()
