import unittest

from source_forward_parity import forward_rows


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


if __name__ == "__main__":
    unittest.main()
