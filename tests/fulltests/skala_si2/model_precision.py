#!/usr/bin/env python3
"""Create an isolated, hash-pinned Float64 diagnostic export, not a new model."""
import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile

import torch


HASHES = {
    "cpu": "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd",
    "cuda": "f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223",
}


def nodes(block):
    for node in list(block.nodes()):
        yield node
        for child in node.blocks():
            yield from nodes(child)


def graphs(model):
    for name, module in model.named_modules():
        for method in module._c._method_names():
            yield name, method, module._c._get_method(method).graph


def inventory(model):
    dtype_args, constants = Counter(), Counter()
    for _, _, graph in graphs(model):
        for node in nodes(graph):
            if node.kind() == "prim::Constant" and node.hasAttribute("value") and node.kindOf("value") == "t":
                constants[str(node.t("value").dtype)] += 1
            if not node.kind().startswith("aten::"):
                continue
            schema = torch._C.parse_schema(node.schema())
            for index, arg in enumerate(schema.arguments):
                if arg.name == "dtype":
                    value = node.inputsAt(index).toIValue()
                    if value is None and node.inputsAt(index).node().kind() != "prim::Constant":
                        raise ValueError("Unresolved dynamic dtype")
                    dtype_args[f"{node.kind()} {value}"] += 1
    return {"parameters": dict(Counter(str(t.dtype) for t in model.parameters())),
            "buffers": dict(Counter(str(t.dtype) for t in model.buffers())),
            "tensor_constants": dict(constants), "dtype_arguments": dict(dtype_args)}


def promote(model):
    # TorchScript encodes dtype as an integer which can share a constant with
    # a shape or index. Replace the dtype operand, never the constant itself.
    original = {name: value.detach().clone() for name, value in model.state_dict().items()}
    model.double()
    for name, value in model.state_dict().items():
        prior = original[name]
        expected = prior.double() if prior.is_floating_point() else prior
        if value.dtype != expected.dtype or not torch.equal(value, expected):
            raise ValueError(f"Parameter/buffer changed beyond exact promotion: {name}")
    changes = []
    for name, method, graph in graphs(model):
        for node in nodes(graph):
            if not node.kind().startswith("aten::"):
                continue
            schema = torch._C.parse_schema(node.schema())
            for index, arg in enumerate(schema.arguments):
                if arg.name != "dtype" or node.inputsAt(index).toIValue() != 6:
                    continue
                if node.kind() not in ("aten::to", "aten::ones"):
                    raise ValueError(f"Unexpected Float32 factory: {node.kind()}")
                replacement = graph.insertConstant(7)
                replacement.node().moveBefore(node)
                node.replaceInput(index, replacement)
                changes.append({"module": name, "method": method, "operation": node.kind(),
                                "argument": arg.name})
        torch._C._jit_pass_lint(graph)
    return changes


def derivative_failures(errors, features):
    failures = []
    for key in features:
        rows = [row for row in errors if row["feature"] == key]
        if [row["step"] for row in rows] != [1e-4, 1e-5, 1e-6] or any(
                not math.isfinite(row[field]) for row in rows
                for field in ("analytic", "finite_difference", "absolute_error")):
            failures.append(key)
            continue
        values = [row["absolute_error"] for row in rows]
        # Allow coarse-step truncation, but require both fine-step bounds and
        # a substantial decrease if the coarse error exceeds those bounds.
        if max(values[1:]) > 1e-8 or (values[0] > 1e-8 and values[1] > values[0]/20):
            failures.append(key)
    return failures


def publish_diagnostic(staged_model, report, output):
    """Publish validated files without replacing any existing model or report."""
    if not report["validation_passed"]:
        raise ValueError("Cannot publish an unvalidated precision model")
    if hashlib.sha256(staged_model.read_bytes()).hexdigest() != report["output_sha256"]:
        raise ValueError("Staged model hash differs from the validation report")
    manifest_text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    os.link(staged_model, output)
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".json",
                                         dir=staged_model.parent) as manifest:
            manifest.write(manifest_text)
            manifest.flush()
            os.link(manifest.name, output.with_suffix(".json"))
    except OSError:
        output.unlink()
        raise


def features(device):
    p = torch.arange(16, dtype=torch.float64, device=device)
    theta = 2 * torch.pi * (p % 8) / 8
    radius = 0.35 + 0.03 * ((p + 1) % 3)
    atoms = torch.tensor([[0., 0., 0.], [0., 0., 1.4]], dtype=torch.float64, device=device)
    coords = atoms[(p // 8).long()] + torch.stack(
        [radius * theta.cos(), radius * theta.sin(), 0.02 * (p % 8 - 3)], dim=1)
    density = (0.25 * (-radius).exp()).repeat(2, 1)
    grad = torch.stack([-0.12 * theta.cos(), -0.12 * theta.sin(), p * 0], dim=0).repeat(2, 1, 1)
    return {"density": density, "grad": grad, "kin": 0.18 * density,
            "grid_coords": coords, "grid_weights": torch.full((16,), 0.04, dtype=torch.float64, device=device),
            "coarse_0_atomic_coords": atoms,
            "atomic_grid_weights": torch.full((16,), 0.05, dtype=torch.float64, device=device),
            "atomic_grid_sizes": torch.tensor([8, 8], dtype=torch.int64, device=device),
            "atomic_grid_size_bound_shape": torch.empty((8, 0), dtype=torch.float64, device=device)}


def evaluate(model, data, derivatives=False):
    inputs = {k: v.detach().clone().requires_grad_(derivatives and v.is_floating_point()
              and k != "atomic_grid_size_bound_shape") for k, v in data.items()}
    energy = (model.get_exc_density(inputs) * inputs["grid_weights"]).sum()
    gradients = {}
    if derivatives:
        keys = [k for k, v in inputs.items() if v.requires_grad]
        values = torch.autograd.grad(energy, [inputs[k] for k in keys], allow_unused=True)
        gradients = {k: v for k, v in zip(keys, values) if v is not None}
    return energy.item(), gradients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--device", choices=HASHES, required=True)
    args = parser.parse_args()
    if args.output.suffix != ".fun":
        parser.error("The diagnostic output must have the .fun suffix")
    if os.path.lexists(args.output) or os.path.lexists(args.output.with_suffix(".json")):
        raise ValueError("Diagnostic output must not already exist")
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    if digest != HASHES[args.device]:
        raise ValueError("Only the pinned upstream Skala-1.1-rev1 export is supported")
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    data = features(args.device)
    with zipfile.ZipFile(args.input) as archive:
        metadata = {name.split("/extra/", 1)[1]: archive.read(name)
                    for name in archive.namelist() if "/extra/" in name}
    if metadata.get("protocol_version") != b"2":
        raise ValueError("Missing Skala protocol metadata")
    model = torch.jit.load(str(args.input), map_location=args.device).eval()
    before = inventory(model)
    baseline_energy, baseline_grad = evaluate(model, data, True)
    # Reload before editing graphs, so no compiled execution plan can be reused.
    model = torch.jit.load(str(args.input), map_location=args.device).eval()
    changes = promote(model)
    after = inventory(model)
    if any(key.endswith(" 6") for key in after["dtype_arguments"]):
        raise ValueError("Float32 dtype argument remains")
    if after["parameters"] != {"torch.float64": 80} or after["buffers"] != {"torch.float64": 43, "torch.int64": 13}:
        raise ValueError("Unexpected parameter/buffer inventory")
    if before["tensor_constants"] != after["tensor_constants"]:
        raise ValueError("Tensor constants unexpectedly changed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata["precision_diagnostic"] = json.dumps({"input_sha256": digest,
        "arithmetic": "float64", "weights": "exact promotion of published float32 weights",
        "purpose": "precision diagnosis only"}).encode()
    with tempfile.TemporaryDirectory(prefix="skala-precision-", dir=args.output.parent) as directory:
        staged_model = Path(directory) / "model.fun"
        torch.jit.save(model, str(staged_model), _extra_files=metadata)
        report = verify_export(staged_model, metadata, after, args.device, data,
                               baseline_energy, baseline_grad, digest, before, changes)
        publish_diagnostic(staged_model, report, args.output)
    print(json.dumps({k: v for k,v in report.items() if k not in ("rewrites", "finite_differences")}, indent=2))
    print("changed dtype operands", len(changes), "max FD error",
          max(r["absolute_error"] for r in report["finite_differences"]))


def verify_export(path, metadata, after, device, data, baseline_energy, baseline_grad,
                  digest, before, changes):
    reloaded_metadata = {key: b"" for key in metadata}
    reloaded = torch.jit.load(str(path), map_location=device,
                              _extra_files=reloaded_metadata).eval()
    if reloaded_metadata != metadata:
        raise ValueError("Serialization lost Skala protocol or diagnostic metadata")
    if inventory(reloaded) != after:
        raise ValueError("Serialization changed the precision inventory")
    energy, gradients = evaluate(reloaded, data, True)
    if not all(math.isfinite(value) for value in (energy, baseline_energy)):
        raise ValueError("Nonfinite diagnostic energies")
    if set(gradients) != set(baseline_grad) or any(
            not torch.isfinite(value).all() for value in (*gradients.values(), *baseline_grad.values())):
        raise ValueError("Missing or nonfinite diagnostic gradients")
    probes = {"density": (1, 2), "grad": (0, 1, 3), "kin": (1, 4),
              "grid_coords": (5, 2), "grid_weights": (6,),
              "coarse_0_atomic_coords": (1, 2), "atomic_grid_weights": (8,)}
    errors = []
    for key, index in probes.items():
        for h in (1e-4, 1e-5, 1e-6):
            minus = {k: v.clone() for k, v in data.items()}
            plus = {k: v.clone() for k, v in data.items()}
            minus[key][index] -= h
            plus[key][index] += h
            fd = (evaluate(reloaded, plus)[0] - evaluate(reloaded, minus)[0]) / (2*h)
            analytic = gradients[key][index].item()
            error = abs(fd-analytic)
            errors.append({"feature": key, "step": h, "analytic": analytic, "finite_difference": fd,
                           "absolute_error": error})
    failures = derivative_failures(errors, probes)
    report = {"scope": "Float64 diagnostic arithmetic with exactly promoted original weights, not a retrained model",
              "input_sha256": digest, "output_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "device": device, "torch": torch.__version__, "before": before, "after": after,
              "rewrites": changes, "baseline_energy": baseline_energy, "diagnostic_energy": energy,
              "energy_difference": abs(energy-baseline_energy),
              "gradient_differences": {k: (v-baseline_grad[k]).abs().max().item() for k,v in gradients.items()},
              "finite_differences": errors, "validation_passed": not failures,
              "failed_features": failures}
    if failures:
        raise ValueError(f"Float64 model derivative failed for {failures}")
    return report


if __name__ == "__main__":
    main()
