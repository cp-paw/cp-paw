import copy
import struct
import unittest

from displace_restart import write_records
import force_mode_parity as modes
import stationarity
from test_force_stationary_fd import protocol


def electronic(step=1):
    return (protocol(step).split("SKALA TOTAL FORCE DIAGNOSTIC")[0]
            + "SKALA TOTAL ENERGY DIAGNOSTIC\nTOTAL ENERGY -4.6938882201616096E+01\n"
            + "NUCLEAR FORCES NOT CALCULATED\nPROGRAM FINISHED\n")


class ForceModeTest(unittest.TestCase):
    def test_energy_only_is_not_a_force(self):
        rows = modes.energy_records(electronic() + electronic(2))
        self.assertEqual([row["step"] for row in rows], [1, 2])
        self.assertNotIn("forces", rows[0])
        self.assertEqual(rows[0]["energy"], -46.938882201616096)

    def test_energy_parser_rejects_incomplete_duplicate_and_nonfinite(self):
        text = electronic()
        bad = [protocol(), text.replace("NUCLEAR FORCES NOT CALCULATED", ""),
               text.replace("TOTAL ENERGY -", "TOTAL ENERGY nan #"),
               text.replace("TOTAL ENERGY -", "TOTAL ENERGY 0\nTOTAL ENERGY -"),
               text.replace("TOTAL ENERGY -", "ATOM 1 0 0 0\nTOTAL ENERGY -"),
               text + "SKALA TOTAL ENERGY DIAGNOSTIC\n",
               text.replace("SKALA TOTAL ENERGY DIAGNOSTIC\n", "")]
        for value in bad:
            with self.subTest(value=value), self.assertRaises(ValueError):
                modes.energy_records(value)

    def test_wave_record_extent(self):
        header = struct.pack("<i", 2) + b"WAVES".ljust(32, b" ")
        data = write_records([b"other", header, b"one", b"two", b"tail"])
        self.assertEqual(modes.wave_records(data), [b"one", b"two"])
        for records in ([b"other"], [header, b"one"], [header, header, b"one", b"two"]):
            with self.assertRaises(ValueError):
                modes.wave_records(write_records(records))

    def test_comparison_checks_electronic_diagnostics_and_layout(self):
        text = electronic()
        data = {"energies": modes.energy_records(text), "trace": stationarity.records(text),
                "bands": stationarity.band_records(text), "scalars": {"MODEL XC ENERGY": [-1.]}}
        self.assertEqual(max(modes.compare(data, data, 1e-10).values()), 0.)
        changes = [lambda d: d["energies"][0].update(energy=0.),
                   lambda d: d["trace"][0].update(maximum=1.),
                   lambda d: d["bands"][0]["bands"][0].update(expectation=1.),
                   lambda d: d["bands"][0]["bands"][0].update(occupation=0.),
                   lambda d: d["bands"].clear(),
                   lambda d: d["scalars"].clear(),
                   lambda d: d["trace"][0].update(step=2)]
        for change in changes:
            other = copy.deepcopy(data)
            change(other)
            with self.assertRaises(ValueError):
                modes.compare(data, other, 1e-10)

    def test_wave_numeric_comparison_preserves_metadata_and_checks_every_payload(self):
        header = struct.pack("<2i9di", 1, 1, *([0.] * 9), 1)
        psi = struct.pack("<8s4i", b"PSI     ", 1, 1, 1, 0)
        grid = struct.pack("<3d3i", 0., 0., 0., 0, 0, 0)
        lam = struct.pack("<8s4i", b"LAMBDA  ", 1, 1, 1, 1)
        values = [header, psi, grid, struct.pack("<2d", 1., 0.), lam, struct.pack("<2d", 2., 0.)]
        self.assertEqual(modes.wave_difference(values, values), 0.)
        changed = list(values)
        changed[-1] = struct.pack("<2d", 2.1, 0.)
        self.assertAlmostEqual(modes.wave_difference(values, changed), 0.1)
        for index, value in ((0, b"bad"), (2, b"bad"), (3, struct.pack("<2d", float("nan"), 0.))):
            changed = list(values)
            changed[index] = value
            with self.assertRaises(ValueError):
                modes.wave_difference(values, changed)

    def test_trajectory_comparison_rejects_changed_metadata(self):
        a = write_records([struct.pack("<idi2d", 1, 0.1, 2, 1., 2.)])
        b = write_records([struct.pack("<idi2d", 1, 0.1, 2, 1.1, 2.)])
        self.assertAlmostEqual(modes.trajectory_difference(a, b), 0.1)
        for record in (struct.pack("<idi2d", 2, 0.1, 2, 1., 2.),
                       struct.pack("<idi2d", 1, 0.1, 2, float("nan"), 2.), b"bad"):
            with self.assertRaises(ValueError):
                modes.trajectory_difference(a, write_records([record]))


if __name__ == "__main__":
    unittest.main()
