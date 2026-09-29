# Periodic Skala integration probes

The CO2, NH3 and urea cells are copied without coordinate changes from
`native-grid-skala-benchmarks/benchmarks/X23-mini/reference`, commit
`ed13cc452542e161a08eb698221b52a56f6c3c47` (the CP2K native-grid Skala
manuscript benchmark repository). Extended XYZ coordinates/cell are in angstrom.

SHA-256:

```
bd245502720bd53126418e5d4ebd2ce520d26bc1cac3df2da66639815bd15717 CO2-solid.xyz
a1b09170f4edac42b77afb568da83ce72636b25d9f6de64c3bd85c77977ef83d NH3-solid.xyz
af4021bb42796f91190c8023ef25ebd533aa77cbc9762542430346271cfbf9bb urea-solid.xyz
```

Run with an opt-in Skala-enabled executable and an external Skala 1.1 model:

```sh
python3 tests/fulltests/skala_crystals/run.py \
  --executable /absolute/path/paw.x --model /absolute/path/model.fun \
  --output /new/directory/crystal-probes --device CUDA
```

The driver performs 180 PBE preparation steps, preserves that restart, and
evaluates Skala with one near-zero-time step. The step count is not an SCF
convergence criterion. `--skala-steps` can request subsequent relaxation.
`--prepare-only` writes inspectable inputs without running them. Existing output
directories are never overwritten. Per-case settings and geometry hashes are
recorded. A nonzero exit or missing completion/adjoint diagnostic is a failure.

This is an integration test, **not** a reference cohesive-energy calculation.
The initial Gamma-only k mesh, 40 Ry cutoff and quadrature must be converged
independently. Reduced `--radial`/`--angular` settings are useful for execution
checks only. The electron-count error is explicitly recorded, not renormalized.
The default does not assert stationary forces, stress or agreement with CP2K.
PAW frozen-core and CP2K GAPW/all-electron energies are not interchangeable.
No D3 correction is included. AlN is not Mani's original case and is not
substituted silently for one of these three crystals.

`--image-shells` controls the compact periodic partition/descriptor support
independently of the quadrature resolution. Results include the constant-field
volume error as well as electron number; neither is normalized away.
`--orientations` adds rotated copies of the selected Lebedev rule and averages
their weights. It does not increase that rule's polynomial exactness.
Electron counts distinguish the reciprocal-space norm plus the PAW overlap
correction, the native pseudo-density grid, and the reconstructed Skala grid.
The first two must agree independently of the atom-grid quadrature. A correct
overlap trace does not certify convergence of the latter; converge cutoff,
radial and angular resolution separately.
See [VALIDATION.md](VALIDATION.md) for the passing execution checks and the
remaining accuracy failures. A successful smoke test is not an EOS validation.
