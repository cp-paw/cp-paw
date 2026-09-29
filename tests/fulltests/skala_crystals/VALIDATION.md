# Joint-source reconstruction checkpoint (2026-09-29)

This is a correctness investigation, not a performance report or a claim that
the AlN equation-of-state issue is fixed. No reference energies were changed.

## Machines and scope

- Spark: GB10, NVHPC 26.5, combined cuBLAS/cuSOLVER/cuFFT build and CUDA Skala.
- Terok: GNU Fortran, OpenBLAS 0.3.33 and CPU-only PyTorch 2.6; no CUDA runtime
  dependency in the CPU executable. Its private static BLAS avoids the two
  BLAS symbol collisions described in the bridge README.
- Isolated remote checkout: `/home/kuehne88/cp-paw-skala-reconstruction`.
  Other CP-PAW checkouts were not modified.
- Reconstruction source checkpoint: `cf051ff`; the subsequent image-shell
  control, volume diagnostic and CPU BLAS fix are part of this follow-up.

## Passing checks

| Check | Result | Scope |
| --- | ---: | --- |
| Joint source/periodic image/partition adjoints, Spark | max error 9.442e-10 | Kernel finite differences |
| Same kernels, Terok GNU | max error 9.431e-10 | Independent compiler |
| Legacy primitive/interpolation kernels | max error 5.930e-10 | Kernel regression |
| GNU CPU-only Si2, 96/17/1 and 200/53/1 | both pass | No CUDA libraries; no preload; private OpenBLAS |
| Ordinary PBE Si2, 180 steps | -7.3657485 H | Existing TYPE=10 reference; passes both without Torch and with private-BLAS CPU Skala build |
| Si2 CPU vs GPU, 200/53/1, 8 k points | total-energy difference 2.727e-7 H | Cold one-step integration; float32 model, not stationary physics |
| Si2 MPI 1 vs 2, 8 k points | energy difference 9.948e-14 H | Same restart and root CUDA inference |
| Same MPI check | force difference 1.086e-8 H/bohr; stress derivative difference 1.570e-8 H | Parallel parity, not force/stress finite-difference validation |

The CPU/GPU energy comparison above used NVHPC/NVPL CPU on Spark. The real
GNU/CPU-only Torch run on Terok is a separate dependency/integration test,
not an interchangeable energy reference for independently initialized orbitals.

## Three molecular crystals

The unchanged structures are from the CP2K native-grid manuscript's X23-mini
data. All use Gamma sampling, 40 Ry, 180 PBE preparation steps and one Skala
snapshot with a near-zero time step. No D3 is included. These settings are
not independently converged, and the PBE preparation count is not a stationarity
criterion. PAW and GAPW energy zeros must not be compared directly.

| Crystal | Radial / angular | Model XC energy (H) | Reconstructed electrons | Expected | Largest operator-contraction error |
| --- | --- | ---: | ---: | ---: | ---: |
| CO2 | 96 / 17 | -89.9758137328103 | 87.9686036171061 | 88 | 2.65e-12 |
| NH3 | 96 / 17 | -32.0499283671351 | 40.0101970994945 | 40 | 3.97e-13 |
| Urea | 96 / 17 | -58.6153041699626 | 64.0055654241851 | 64 | 6.27e-13 |
| CO2 | 200 / 53 | -89.9833667782936 | 88.0076516249965 | 88 | 2.59e-11 |

Every row completed successfully on Spark. The CO2 XC energy changes by about
0.00755 H upon refinement: these are deliberately not accepted reference
energies. Logs/results are in `validation/crystals-integration` and
`validation/co2-grid-refinement` on Spark.

A repeat with the volume diagnostic and exact-zero pair-loop shortcuts gives
the same three coarse-grid XC energies to all digits printed above. All three
operator checks pass. Its constant-field diagnostics are:

| Crystal, 96/17/1 | Relative volume error |
| --- | ---: |
| CO2 | -0.2523% |
| NH3 | -0.4837% |
| Urea | +0.1199% |

These errors mix angular/radial and finite-image effects; they must not be
attributed to the image layout alone. The runs completed in
`validation/crystals-volume-check`; their provenance files record structure,
executable and model hashes. All tests at this checkpoint have finished.

## Accuracy failure isolated without Skala

The constant-field volume integral on the Si2 primitive cell fails even
without density reconstruction or model inference. Exact volume is
270.011394 bohr^3. With 200 radial points, angular exactness 53 and one
orientation:

| Fixed image shells | Integrated volume (bohr^3) | Relative error |
| --- | ---: | ---: |
| 1 | 276.114991162 | 2.2605% |
| 2 | 270.874984776 | 0.3198% |
| 3 | 270.247379140 | 0.0874% |

Doubling radial points alone barely changes the Si2 particle count
(28.1804442 to 28.1804135); it cannot remove this image-layout error.
The independent `partition_measure.x` probe now exposes this failure. The
production path reports its partition volume and warns on relative error
above 1e-4. `IMAGESHELLS` is an independent convergence parameter; its default
of one shell must not be interpreted as converged. Neither densities nor
weights are rescaled to fit a reference.

The applied Si2 model with two image shells also completes on Spark. Its
electron count is 28.0256553536904, versus 28.1804442 with one shell (exact
count 28). Model XC energy is -41.4502546897069 H and its largest operator
contraction error is 2.7e-15. This isolates a substantial image-layout effect;
it is not a converged energy. The three-shell volume probe still fails the
1e-4 relative-accuracy target.

## Still required

The following follow-up supersedes the *partition implementation* in the
historical tables above; those earlier energies are retained as diagnostic
history, not updated scientific references.

### Compact periodic partition and independent trace

The kernel introduced in `8617578` is now connected to both energy weights
and the separate self-image descriptor window, including atom and cell
derivatives. See the reconstruction unit-test README for the mathematical
definition and its finite-difference tests. It is a revised finite-shell
discretization; energy/descriptor shell convergence is still necessary.

Model-free Si2 constant-field checks:

| Shells / radial / angular exactness / orientations | Relative volume error | First cosine integral (bohr^3) | 1e-4 relative tolerance |
| --- | ---: | ---: | --- |
| 1 / 200 / 53 / 1 | 5.606297506e-5 | -0.006812607460 | Pass |
| 1 / 400 / 53 / 1 | 5.606297863e-5 | -0.006812607749 | Pass |
| 1 / 200 / 53 / 3 | 2.600373546e-5 | -0.003484566173 | Pass |
| 2 / 200 / 53 / 1 | 1.809252458e-4 | -0.021658174220 | **Fail** |

Increasing radial resolution alone is ineffective here. Increasing the image
support does not guarantee a smaller finite-quadrature error at a fixed
angular rule. The failed two-shell test is retained, not used to select a
convenient shell default. More angular refinement and model-energy convergence
remain required.

Applied Skala Si2 checks use the periodic two-atom diamond primitive cell,
eight k points, 200 radial points, angular exactness 53 (974 directions),
one orientation and one image shell. These are cold, near-zero-time-step
integration probes, not stationary electronic states:

| Build / EPWPSI | All-electron overlap trace | Reconstructed grid electrons | Grid minus trace | Model XC energy (H) |
| --- | ---: | ---: | ---: | ---: |
| Spark GPU / 20 | 28.0000000000000 | 28.0005725362595 | 5.725362595e-4 | -41.4437984610168 |
| Terok GNU CPU-only / 20 | 27.9999999906291 | 28.0005725268547 | 5.725362256e-4 | -41.4437987593554 |
| Spark GPU / 40 | 28.0000000000000 | 28.0005676185619 | 5.676185619e-4 | -41.4431806483253 |

The reciprocal pseudo-wavefunction norm plus `Tr(D DeltaS)` gives the
valence count; frozen-core occupations are added separately. The native
pseudo-density grid agrees with the reciprocal norm to 2e-14 on CPU and
5e-14 or better on GPU. The CPU occupation-trace difference is -9.371e-9,
consistent with the legacy iterative orthogonalizer's `max|S-I| < 1e-8`
criterion; the recommended GPU build uses the initial Cholesky path. Tests
therefore bound the occupation-trace difference by that tolerance times the
valence occupation sum, not by a tighter arbitrary absolute threshold.

At EPWPSI=20 the CPU/GPU total-energy difference is 3.238e-7 H, with the
float32 model and independently initialized orbitals. The largest displayed
operator-contraction discrepancy across these three runs is 1.821e-14.
The ordinary 180-step PBE Si2 regression also passes with the new CPU binary
at its unchanged -7.3657485 H reference.

The cutoff change barely changes the electron integration error in this
particular probe. It does change the orbital representation, so the two XC
energies alone are not a fixed-density interpolation or SCF convergence test.
No density renormalization was applied. The native-grid cutoff, atom-grid
radial/angular accuracy and image support must be converged independently.

Remote evidence directories under the isolated checkout:
`validation/si2-periodic-gpu-r2`, `validation/si2-periodic-cpu-r2`,
`validation/si2-periodic-cutoff40`, and `validation/pbe-periodic-regression`.
The first trial in `si2-periodic-{cpu,gpu}` stopped on a setup-selection
guard in the new diagnostic; that programming error was fixed before the
successful `-r2` runs. SHA-256 provenance:

```text
Spark executable cbfe373c3dca7a241c87e5f94e0e6dbca4fb4533de61292609156246745d98d3
Terok executable 6e0070c8d274f3b4cc0fc18576b1224dc923f728a42e081693d4fffcbf656550
CUDA model f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223
CPU model 7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd
paw_skala.f90 29b5335727e088a5aaa25e0c492800c55f9865614d5f2fdc14286884518e478b
paw_skala_partition.f90 e534c00660392f2025c49de6a378ab49eea8e467a53be4caf6ef49e26615d549
paw_waves1.f90 d11057f68d79ee531705a92ff5a03f4131e7f83ee77d6d30ac1f73a00e752b66
```

The new crystal driver records all three electron-count routes; its four
Python tests pass. The three molecular-crystal calculations in the earlier
table have **not yet been repeated with this new partition**. The new
end-to-end MPI check is also pending. The NVHPC/LibTorch executable still
emits a multiple-OpenMP-runtime warning; these checks used one OpenMP thread.

### Remaining acceptance gates

1. Resolve/converge the periodic energy partition and descriptor image window,
   as well as radial/angular and native-grid interpolation errors.
2. Relax the electronic states with the corrected applied functional.
3. Validate full PAW forces and stress against stationary-state total-energy
   finite differences, beyond kernel and MPI checks.
4. Reproduce the AlN equation of state with documented PAW setups, k sampling,
   quadrature, restart provenance and separate dispersion settings. Mani's
   original CP-PAW inputs were not available.
5. Audit the nonzero one-center `TAU-KINETIC DIFFERENCE` diagnostic separately:
   its gradient-form positive tau and the setup kinetic matrix need matching
   radial domains, boundary terms and relativistic conventions before claiming
   equality. The current Si2 differences (0.00786 and 0.00399 H) are recorded,
   not certified as errors or dismissed as harmless.

Until these gates pass, the joint reconstruction remains experimental and is
not ready for a scientific acceptance claim or merge.
