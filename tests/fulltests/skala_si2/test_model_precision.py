import importlib.util
import hashlib
import json
import tempfile
import unittest
from pathlib import Path


HAVE_TORCH = importlib.util.find_spec("torch") is not None
if HAVE_TORCH:
    import torch
    import model_precision

    class SharedDtypeAndShape(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.weight = torch.nn.Parameter(torch.arange(6, dtype=torch.float32))
            self.register_buffer("offset", torch.full((6,), 0.25, dtype=torch.float32))
            self.register_buffer("indices", torch.arange(6, dtype=torch.int64))

        def forward(self, x):
            value = x.to(torch.float32) + self.weight + self.offset
            if x[0] > 0:
                value = value + torch.ones(6, dtype=torch.float32)
            else:
                value = value - torch.ones(6, dtype=torch.float32)
            return value, self.indices

    class DynamicDtype(torch.nn.Module):
        def forward(self, x):
            return torch.ones(6, dtype=x.dtype)


@unittest.skipUnless(HAVE_TORCH, "precision export tests require PyTorch")
class ModelPrecisionTest(unittest.TestCase):
    def test_exact_promotion_preserves_integer_shape_and_nested_blocks(self):
        model = torch.jit.script(SharedDtypeAndShape())
        before = model_precision.inventory(model)
        self.assertIn("aten::to 6", before["dtype_arguments"])
        changes = model_precision.promote(model)
        self.assertGreaterEqual(len(changes), 3)
        after = model_precision.inventory(model)
        self.assertNotIn("aten::to 6", after["dtype_arguments"])
        self.assertNotIn("aten::ones 6", after["dtype_arguments"])
        self.assertEqual(after["parameters"], {"torch.float64": 1})
        self.assertEqual(after["buffers"], {"torch.float64": 1, "torch.int64": 1})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "precision-check.pt"
            torch.jit.save(model, str(path))
            reloaded = torch.jit.load(str(path))
            self.assertEqual(model_precision.inventory(reloaded), after)
            for sign in (-1., 1.):
                x = torch.full((6,), sign * 0.125, dtype=torch.float64)
                value, indices = reloaded(x)
                expected = x + torch.arange(6, dtype=torch.float64) + 0.25 + sign
                self.assertEqual(value.dtype, torch.float64)
                self.assertTrue(torch.equal(value, expected))
                self.assertEqual(indices.dtype, torch.int64)
                self.assertEqual(indices.tolist(), list(range(6)))

    def test_dynamic_dtype_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "dynamic dtype"):
            model_precision.inventory(torch.jit.script(DynamicDtype()))

    def test_derivative_gate_requires_finite_complete_fine_step_evidence(self):
        rows = [{"feature": "kin", "step": h, "analytic": 1.,
                 "finite_difference": 1. + error, "absolute_error": error}
                for h, error in zip((1e-4, 1e-5, 1e-6), (6e-8, 6e-10, 6e-12))]
        self.assertEqual(model_precision.derivative_failures(rows, ["kin"]), [])
        self.assertEqual(model_precision.derivative_failures(rows[:2], ["kin"]), ["kin"])
        for invalid in (float("nan"), float("inf"), 2e-8):
            changed = [dict(row) for row in rows]
            changed[-1]["absolute_error"] = invalid
            self.assertEqual(model_precision.derivative_failures(changed, ["kin"]), ["kin"])
        rows[1]["absolute_error"] = 5e-9
        self.assertEqual(model_precision.derivative_failures(rows, ["kin"]), ["kin"])

    def test_publish_requires_validation_and_never_replaces_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            staged = root / "staged.fun"
            staged.write_bytes(b"diagnostic")
            output = root / "output.fun"
            report = {"output_sha256": hashlib.sha256(staged.read_bytes()).hexdigest(),
                      "validation_passed": False}
            with self.assertRaisesRegex(ValueError, "unvalidated"):
                model_precision.publish_diagnostic(staged, report, output)
            self.assertFalse(output.exists())
            report["validation_passed"] = True
            model_precision.publish_diagnostic(staged, report, output)
            self.assertEqual(output.read_bytes(), staged.read_bytes())
            self.assertEqual(json.loads(output.with_suffix(".json").read_text()), report)
            with self.assertRaises(FileExistsError):
                model_precision.publish_diagnostic(staged, report, output)
            self.assertEqual(output.read_bytes(), b"diagnostic")
            output.unlink()
            output.with_suffix(".json").write_text("preserve existing report")
            with self.assertRaises(FileExistsError):
                model_precision.publish_diagnostic(staged, report, output)
            self.assertFalse(output.exists())
            self.assertEqual(output.with_suffix(".json").read_text(), "preserve existing report")


if __name__ == "__main__":
    unittest.main()
