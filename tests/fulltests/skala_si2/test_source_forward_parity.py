import csv
from pathlib import Path
import tempfile
import unittest

from source_forward_parity import check_resident_profile, check_reuse, forward_rows, reused_rows


class ForwardCoverageTests(unittest.TestCase):
    def test_exact_evaluation_count(self):
        text = "SOURCE GPU FORWARD ROWS 0\nSOURCE GPU FORWARD ROWS 35420\n"
        self.assertEqual(forward_rows(text, 2), [0, 35420])

    def test_missing_or_extra_evaluation(self):
        for text in ("", "SOURCE GPU REVERSE ROWS 20\n",
                     "SOURCE GPU FORWARD ROWS 1\nSOURCE GPU FORWARD ROWS 2\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                forward_rows(text, 1)

    def test_invalid_coverage(self):
        for value in ("-1", "1.5", "NaN", "Infinity", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                forward_rows("SOURCE GPU FORWARD ROWS " + value, 1)

    def test_reused_rows_have_their_own_label(self):
        self.assertEqual(reused_rows("SOURCE GPU REUSED ROWS 0\nSOURCE GPU REUSED ROWS 20\n", 2), [0, 20])
        with self.assertRaises(ValueError):
            reused_rows("SOURCE GPU FORWARD ROWS 20\n", 1)

    def test_residency_modes(self):
        check_reuse([20, 20], [0, 20], "full")
        check_reuse([20, 20], [0, 10], "bounded")
        check_reuse([20, 20], [0, 0], "off")
        for hits, mode in (([1, 20], "full"), ([0, 21], "full"), ([0, 0], "full"),
                           ([0, 20], "bounded"), ([0, 0], "bounded"), ([0, 1], "off")):
            with self.subTest(hits=hits, mode=mode), self.assertRaises(ValueError):
                check_reuse([20, 20], hits, mode)


class ResidencyProfileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.rows = [
            ["ACC_COPY_SKALA_FWD_GEOM", 10, 2, 12, 4, 1, 2060e-9],
            ["ACC_COPY_SKALA_FWD_INPUT", 10, 2, 12, 4, 2, 64e-9],
            ["SKALA_SOURCE_FWD_DEVICE", 10, 2, 12, 4, 2, 0],
            ["ACC_COPY_SKALA_FWD_OUT", 10, 2, 12, 4, 2, 1600e-9],
            ["ACC_PRESENT_SKALA_FWD_GEOM", 10, 2, 12, 4, 2, 0]]

    def write(self):
        with (self.work / "source_forward_profile.csv").open("w") as stream:
            writer = csv.writer(stream)
            writer.writerow(["op", "n1", "n2", "n3", "n4", "calls", "gbyte"])
            writer.writerows(self.rows)

    def test_first_upload_and_warm_hit(self):
        self.write()
        result = check_resident_profile(self.work, 1, [10, 10], [0, 10])
        self.assertEqual(result["ACC_COPY_SKALA_FWD_GEOM"]["rows"], 10)
        self.assertEqual(result["ACC_PRESENT_SKALA_FWD_GEOM"]["rows"], 20)

    def test_missing_electronic_input(self):
        self.rows.pop(1)
        self.write()
        with self.assertRaises(ValueError):
            check_resident_profile(self.work, 1, [10, 10], [0, 10])

    def test_incorrect_payload(self):
        for index in (0, 3, 4):
            before = self.rows[index][-1]
            self.rows[index][-1] += 100e-9
            self.write()
            with self.subTest(index=index), self.assertRaises(ValueError):
                check_resident_profile(self.work, 1, [10, 10], [0, 10])
            self.rows[index][-1] = before

    def test_false_warm_hit(self):
        self.write()
        with self.assertRaises(ValueError):
            check_resident_profile(self.work, 1, [10, 10], [0, 0])

    def test_missing_rank(self):
        self.write()
        with self.assertRaises(ValueError):
            check_resident_profile(self.work, 2, [10, 10], [0, 10])


if __name__ == "__main__":
    unittest.main()
