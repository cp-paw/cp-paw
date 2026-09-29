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

1. Resolve/converge the periodic energy partition and descriptor image window,
   as well as radial/angular and native-grid interpolation errors.
2. Relax the electronic states with the corrected applied functional.
3. Validate full PAW forces and stress against stationary-state total-energy
   finite differences, beyond kernel and MPI checks.
4. Reproduce the AlN equation of state with documented PAW setups, k sampling,
   quadrature, restart provenance and separate dispersion settings. Mani's
   original CP-PAW inputs were not available.

Until these gates pass, the joint reconstruction remains experimental and is
not ready for a scientific acceptance claim or merge.
