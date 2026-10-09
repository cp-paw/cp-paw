"""Regression checks for shared CP-PAW binary-record helpers."""

import struct
import unittest

import paw_records as records


def section(name, count=1):
    return struct.pack("<i32s", count, name.encode("ascii").ljust(32))


class RecordTests(unittest.TestCase):
    def test_roundtrip(self):
        payloads = [b"", b"test", struct.pack("<2d", 1., -2.)]
        self.assertEqual(records.read_records(records.write_records(payloads)), payloads)

    def test_invalid_markers(self):
        data = records.write_records([b"test"])
        for bad in (data[:3], data[:-1], data[:-4] + struct.pack("<I", 3)):
            with self.subTest(data=bad), self.assertRaises(ValueError):
                records.read_records(bad)

    def test_displacement_preserves_other_records(self):
        rows = [section("ATOMS"), struct.pack("<i", 2), b"metadata",
                struct.pack("<6d", *range(6)), struct.pack("<6d", *range(6, 12)),
                section("CELL"), struct.pack("<27d", *range(27)),
                b"untouched electronic payload"]
        before = records.write_records(rows)
        shifted = records.displace_geometry(before, 2, 1, 0.125)
        after = records.read_records(shifted)
        for index in (0, 1, 2, 5, 6, 7):
            self.assertEqual(rows[index], after[index])
        for old, new in zip(records.geometry(before)["positions"],
                            records.geometry(shifted)["positions"]):
            self.assertEqual([0., 0., 0., 0.125, 0., 0.], [b-a for a, b in zip(old, new)])
        translated = records.geometry(records.displace_geometry(
            before, 1, 2, 0.125, rigid_translation=True))
        for old, new in zip(records.geometry(before)["positions"], translated["positions"]):
            self.assertEqual([0., 0.125, 0., 0., 0.125, 0.], [b-a for a, b in zip(old, new)])

    def test_trajectory_validation(self):
        def trajectory(time, value):
            return records.write_records([struct.pack("<idid", 1, time, 1, value)])
        self.assertEqual(records.trajectory_difference(trajectory(1., 2.), trajectory(1., 2.125)),
                         0.125)
        for bad in (trajectory(2., 2.), trajectory(1., float("nan"))):
            with self.subTest(data=bad), self.assertRaises(ValueError):
                records.trajectory_difference(trajectory(1., 2.), bad)

    def test_wave_validation(self):
        waves = [struct.pack("<2i9di", 1, 1, *range(9), 1),
                 struct.pack("<8s4i", b"PSI     ", 1, 1, 1, 0), bytes(36),
                 struct.pack("<2d", 2., 0.),
                 struct.pack("<8s4i", b"LAMBDA  ", 1, 1, 1, 1),
                 struct.pack("<2d", 3., 0.)]
        data = records.write_records([section("WAVES", len(waves)), *waves])
        self.assertEqual(records.wave_records(data), waves)
        changed = list(waves)
        changed[3] = struct.pack("<2d", 2.125, 0.)
        self.assertEqual(records.wave_difference(waves, changed), 0.125)
        changed[3] = struct.pack("<2d", float("nan"), 0.)
        with self.assertRaises(ValueError):
            records.wave_difference(waves, changed)
        with self.assertRaises(ValueError):
            records.wave_records(records.write_records([section("WAVES", 99), *waves]))


if __name__ == "__main__":
    unittest.main()
