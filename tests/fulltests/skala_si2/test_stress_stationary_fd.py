import copy
import json
from pathlib import Path
import struct
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from deform_restart import deform, read_records, write_records
import force_stationary_fd as forces
import stress_stationary_fd as stress
from test_force_stationary_fd import restart
from test_stress_mode_parity import stress_protocol


def wave_restart():
    header = struct.pack("<2i9di", 1, 1, *([0.]*9), 1)
    psi = struct.pack("<8s4i", b"PSI     ", 1, 1, 1, 0)
    grid = struct.pack("<3d3i", 0., 0., 0., 0, 0, 0)
    lam = struct.pack("<8s4i", b"LAMBDA  ", 1, 1, 1, 1)
    values = [header, psi, grid, struct.pack("<2d", 1., 0.), lam, struct.pack("<2d", 2., 0.)]
    return write_records([struct.pack("<i", len(values))+b"WAVES".ljust(32, b" "), *values])


class StationaryStressTest(unittest.TestCase):
    def test_affine_transform_changes_all_geometry_levels_only(self):
        data = restart()+wave_restart()
        matrix = stress.direction("xy")
        result = deform(data, [0.02*x for row in matrix for x in row])
        before, after = read_records(data), read_records(result)
        for i, (old, new) in enumerate(zip(before, after)):
            if i == 1:
                self.assertEqual(struct.unpack("<27d", new),
                    (1., .01, 0., .01, 1., 0., 0., 0., 1.)*3)
            elif i in (5, 6):
                self.assertEqual(struct.unpack("<6d", new), (0., 0., 0., 2.03, 3.02, 4.))
            else:
                self.assertEqual(old, new)
        self.assertEqual(deform(data, [0.]*9), data)

    def test_bad_strains_and_malformed_geometry_fail(self):
        for strain in ([0.]*8, [float("nan")]+[0.]*8, [-1.]+[0.]*8, [-2.]+[0.]*8,
                       [1e308, 1e308]+[0.]*7):
            with self.subTest(strain=strain), self.assertRaises(ValueError):
                deform(restart(), strain)
        records = read_records(restart())
        for data in (b"", restart()[:-1], write_records(records+records),
                     write_records(records[:3]), write_records([records[0], b"bad", *records[2:]])):
            with self.assertRaises(ValueError):
                deform(data, [0.]*9)

    def test_basis_fingerprint_ignores_cell_and_wave_values_but_not_g_vectors(self):
        data = wave_restart()
        original = stress.basis_signature(data)
        records = read_records(data)
        struct.pack_into("<d", records[1], 8, 1.)
        struct.pack_into("<d", records[4], 0, .99)
        self.assertEqual(stress.basis_signature(write_records(records)), original)
        struct.pack_into("<i", records[3], 24, 1)
        self.assertNotEqual(stress.basis_signature(write_records(records)), original)
        struct.pack_into("<d", records[4], 0, float("nan"))
        with self.assertRaises(ValueError):
            stress.basis_signature(write_records(records))

    def test_tensor_contraction_sign_and_shear_convention(self):
        center = forces.analyze(stress_protocol(), 2, 1, 1, 1e-7, 1e-7)
        center.update(final_stress=stress.stress_records(stress_protocol())[-1],
                      basis_signature="same", grid_signature=[570, 4096])
        for name, value in (("isotropic", 6.), ("xx", 1.), ("xy", .1), ("xz", .2), ("yz", .3)):
            minus, plus = json.loads(json.dumps(center)), copy.deepcopy(center)
            minus["final"]["energy"] -= value*.001
            plus["final"]["energy"] += value*.001
            row = stress.compare(center, minus, plus, .001, name, 1e-10)
            self.assertTrue(row["passed"])
            self.assertAlmostEqual(row["analytic_hartree"], value)
            self.assertIsNone(stress.compare(center, minus, plus, .001, name)["passed"])
        with self.assertRaises(ValueError):
            stress.direction("bad")

    def test_nonstationary_changed_basis_and_occupations_block_comparison(self):
        center = forces.analyze(stress_protocol(), 2, 1, 1, 1e-7, 1e-7)
        center.update(final_stress=stress.stress_records(stress_protocol())[-1],
                      basis_signature="same", grid_signature=[570, 4096])
        for change in (dict(stationary=False), dict(occupations=[]), dict(basis_signature="new"),
                       dict(grid_signature=[571, 4096])):
            plus = dict(center, **change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                stress.compare(center, center, plus, .001, "xx")

    def test_relaxation_retains_failed_blocks_and_requires_complete_stress(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = {key: root/key for key in ("restart", "model", "structure", "executable")}
            for path in inputs.values():
                path.write_bytes(b"fixture")
            inputs["restart"].write_bytes(restart()+wave_restart())
            args = SimpleNamespace(output=root, atom=1, axis=1, dt=5, block_steps=1,
                max_blocks=2, last=1, residual_tolerance=1e-7, commutator_tolerance=1e-7,
                device="CPU", radial_points=96, lebedev_exactness=17, cutoff=20.,
                mass=25., mass_g2=.3166286988823056, friction=.1, timeout=1, stress=True)
            def run(command, cwd, **kwargs):
                text = stress_protocol()
                if cwd.name == "block-1":
                    text = text.replace("1e-8", "1e-2")
                (cwd/"si2.prot").write_text(text+"#(G-VECTORS FOR DENSITY)...: 570\nGRID POINTS 4096\n")
            with patch.object(forces.subprocess, "run", side_effect=run):
                _, result = stress.run_leg("xx", [.001]+[0.]*8, args, inputs, {})
            self.assertTrue(result["stationary"])
            self.assertEqual([x["stationary"] for x in result["blocks"]], [False, True])
            self.assertEqual(result["final_stress"]["tensor"][0][0], 1.)
            control = (root/"xx/block-1/si2.cntl").read_text()
            self.assertIn("STRESS=T", control)
            self.assertIn("!CELL MOVE=F", control)
            from test_force_mode_parity import electronic
            args.electronic_warmup_steps = 1
            def warmup(command, cwd, **kwargs):
                if cwd.name.startswith("electronic"):
                    (cwd/"si2.prot").write_text(electronic()+
                        "#(G-VECTORS FOR DENSITY)...: 570\nGRID POINTS 4096\n")
                else:
                    run(command, cwd, **kwargs)
            with patch.object(forces.subprocess, "run", side_effect=warmup):
                _, result = stress.run_leg("warmup", [0.]*9, args, inputs, {})
            self.assertTrue(result["stationary"])
            self.assertEqual(len(result["electronic_blocks"]), 2)
            self.assertEqual(result["final_stress"]["tensor"][0][0], 1.)
            def changed_grid(command, cwd, **kwargs):
                warmup(command, cwd, **kwargs)
                if cwd.name.startswith("electronic"):
                    path = cwd/"si2.prot"
                    path.write_text(path.read_text().replace("GRID POINTS 4096", "GRID POINTS 8192"))
            with patch.object(forces.subprocess, "run", side_effect=changed_grid):
                _, result = stress.run_leg("grid-changed", [0.]*9, args, inputs, {})
            self.assertFalse(result["stationary"])
            self.assertIn("Grid size changed", result["failure"])
            args.electronic_warmup_steps = 0
            def missing(command, cwd, **kwargs):
                run(command, cwd, **kwargs)
                file = cwd/"si2.prot"
                file.write_text(file.read_text().replace("TOTAL D E / D STRAIN 0.2 0.3 3.0\n", ""))
            with patch.object(forces.subprocess, "run", side_effect=missing):
                _, result = stress.run_leg("missing", [0.]*9, args, inputs, {})
            self.assertFalse(result["stationary"])
            self.assertIn("stress", result["failure"].lower())


if __name__ == "__main__":
    unittest.main()
