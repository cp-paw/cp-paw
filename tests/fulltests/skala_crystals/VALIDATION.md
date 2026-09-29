# Joint-source reconstruction checkpoint (2026-09-29)

This is a correctness investigation, not a performance report or a claim that
the AlN equation-of-state issue is fixed. No reference energies were changed.

## Machines and scope

- Spark: 20 CPU cores and GB10, NVHPC 26.5, combined cuBLAS/cuSOLVER/cuFFT
  build with CUDA Skala, plus CPU-mode validation with one/eight MPI ranks.
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
| 1 / 200 / 59 / 1 | 2.312483851e-5 | -0.002819803756 | Pass |
| 1 / 200 / 65 / 1 | 9.684943826e-6 | -0.001185767227 | Pass |
| 2 / 200 / 65 / 1 | 4.166476066e-5 | -0.005028693319 | Pass |

Increasing radial resolution alone is ineffective here. Increasing the image
support does not guarantee a smaller finite-quadrature error at a fixed
angular rule. The failed two-shell test is retained, not used to select a
convenient shell default. More angular refinement and model-energy convergence
remain required. The extended order-65 rule reduces the two-shell integration
error below the stated tolerance; it does not remove the need for separate
image-shell and energy/force convergence tests.

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

At this checkpoint the crystal repeats and current-build MPI checks were
still pending. The completed follow-up checks are recorded below. The
NVHPC/LibTorch executable still emits a multiple-OpenMP-runtime warning;
these checks used one OpenMP thread without suppressing the warning.

### Molecular crystals with the current partition

All three structures have now completed applied Skala snapshots with the
current periodic partition and electronic diagnostics. CO2 and NH3 use the
GNU CPU-only build on Terok; urea uses Spark's NVHPC/NVPL binary in explicit
CPU mode (`CPPAW_GPU_MODE=off`, `DEVICE='CPU'`, no visible CUDA devices).
All use eight MPI ranks, one OpenMP/OpenBLAS thread per rank, the same CPU
model export, 96 radial points, angular order 17, one orientation and one
image shell. The unchanged preparation is 180 PBE steps, Gamma sampling,
40 Ry, no D3, followed by one near-zero-time-step Skala snapshot.

| Crystal / CPU host | Model XC energy (H) | All-electron trace | Grid minus trace | Relative volume error | Occupied maximum residual (H) |
| --- | ---: | ---: | ---: | ---: | ---: |
| CO2 / Terok | -89.9746595603607 | 88 | -3.879892690e-2 | -0.45066585% | 8.520398811e-2 |
| NH3 / Terok | -32.0491171650467 | 40 | +5.914619114e-3 | -0.64071940% | 7.563157745e-2 |
| Urea / Spark | -58.6138316838099 | 64 | -1.985967666e-3 | -0.09182597% | 7.786199297e-2 |

The native pseudo-density grid and reciprocal norm agree within 2.843e-14
electrons. The independent PAW valence trace agrees with occupations within
1.422e-14 electrons. The largest tested adjoint/operator-contraction error is
2.843e-14, overlap error 1.111e-15 and Hamiltonian Hermiticity error
1.321e-15 H. All three completed protocols pass diagnostic-consistency checks.
None passes the 1e-4 relative-volume target at this coarse quadrature; their
Skala electronic states are also not stationary. These results therefore
remain integration probes, not physical reference energies or an AlN EOS fix.
No weights, densities or energies were adjusted to remove the errors.

Evidence: Terok `validation/crystals-periodic-cpu-mpi8/{CO2,NH3}` and Spark
`validation/crystal-urea-spark-cpu-mpi8/urea`. Their executable/model hashes
are the same as the CPU MPI checks below; per-case `provenance.json` files
also record structure hashes and settings. Spark's urea protocol SHA-256 is
`8871324e691b477ace386c8938fab1c3290c6381a129876a2558c061dd01e516`.

The original serial Terok CO2 attempt was intentionally stopped to switch to
MPI and is not counted as a pass. After CO2 and NH3 completed, the redundant
Terok urea attempt was stopped because urea was already running on Spark;
the Terok batch exit 143 is a cancellation, not three completed cases.
All partial logs are retained. Spark's numerical urea run finished normally,
but its initially loaded Python harness exited 1 on a field-name collision
between `PARTITION VOLUME` and `PARTITION VOLUME RELATIVE ERROR`. The parser
was corrected, and all 22 scalar fields in each completed crystal protocol
were rechecked without rerunning or modifying the calculations. Its eight
unit tests now cover both that collision and hidden/nonfinite earlier values.
The successful Spark reanalysis is saved separately in
`validation/crystal-urea-spark-cpu-mpi8-recheck.log`; the initial failure log
has not been overwritten.

### Higher Lebedev rules

Orders 59 (1202 points) and 65 (1454 points) were added without changing the
default 53. A requested minimum of 64 selects order 65. The independent
polynomial tests pass with both GNU/Terok and NVHPC/Spark; see
[`LEBEDEV.md`](../../unittests/skala_reconstruction/LEBEDEV.md) for coefficient
provenance, moment errors and the model-free angular sweep.

The applied Si2 check was repeated on Spark with requested exactness 64,
200 radial points, one orientation, one image shell, EPWPSI=20 and eight
k points. The protocol confirms actual order 65 and 1454 directions. Like
the order-53 comparison, this is a cold near-zero-time-step integration probe,
not an electronically stationary state:

| Actual angular order | Grid minus all-electron trace | Relative volume error | Model XC energy (H) |
| ---: | ---: | ---: | ---: |
| 53 | 5.725362595e-4 | 5.606297506e-5 | -41.4437984610168 |
| 65 | 8.741174813e-5 | 9.684943808e-6 | -41.4435426565372 |

The order-65 overlap trace is 28 electrons; the integrated atom-grid count
is 28.0000874117481. The latter error is about 6.55 times smaller than at
order 53, without density normalization. The model XC energy changes by
2.558044796e-4 H, so these two grids alone do not establish energy convergence.
The maximum reported adjoint/operator-contraction error is 4.781e-15, all
tested strain tensors are finite, and the integration test passes. These
checks are not stationary-state force/stress finite differences.

Evidence: Spark `validation/si2-periodic-lebedev64`, with executable SHA-256
`ef88b93117beed97c3d36e159aba812d52670beb70f0e4623a059b881b44e153` and the
same CUDA model as above. The default remains 53 pending a broader convergence
study; the new orders are selectable on both CPU and GPU.

### Electronic stationarity and CPU MPI

The new `CHECK=T` diagnostic evaluates `H Psi - S Psi (Psi^dagger H Psi)`,
the occupied RMS and maximum norm, the occupation commutator, overlap error,
and Hamiltonian Hermiticity. It is read-only and follows the Hamiltonian
selected by `APPLY`. Passing the residual alone does not establish correct
occupations, orthonormality, a global minimum, or quadrature convergence.

A new eight-k-point PBE warm start (300 steps, 20 Ry, `CDUAL=2`) matches this
Si2 structure. The subsequent diagnostic runs use `DT=1e-6`, 96 radial points,
order 17, one orientation and one shell on Terok GNU/CPU-only. These deliberately
coarse quadratures test instrumentation, not physical convergence:

| State / Hamiltonian | Occupied RMS (H) | Occupied maximum (H) | Occupation commutator (H) | Overlap error |
| --- | ---: | ---: | ---: | ---: |
| Cold / Skala | 1.729613374 | 1.889739732 | 0.3129584933 | 2.001e-9 |
| PBE warm restart / PBE (`APPLY=F`) | 2.123543829e-8 | 6.449052927e-8 | 1.448571787e-7 | 2.442e-15 |
| Same PBE restart / Skala (`APPLY=T`) | 9.426973264e-2 | 1.405628165e-1 | 1.139137897e-2 | 2.442e-15 |

Hamiltonian Hermiticity errors are below 6e-16 H. The independent PBE probe
passes explicit 1e-6 H maximum-residual and commutator limits. The switched
Skala probe fails those limits, as it should. No energy-change or step-count
criterion was substituted for electronic stationarity. The parser has nine
tests, including missing/nonfinite diagnostics, earlier-step failures, a
stationary subspace with wrong occupations, and final-window handling.

The updated CPU-only MPI executable also passes 1/2-rank parity from the same
restart: energy difference 2.295e-12 H, force difference 3.665e-10 H/bohr,
strain-derivative difference 9.001e-10 H, smooth-operator norm difference
1.135e-8, and maximum difference among the electronic diagnostics 1.740e-10.
This closes the CPU MPI consistency check. The GPU MPI check is reported
separately below; stationary-state force/stress finite differences remain open.

Evidence: Terok `validation/si2-stationarity-{cold,pbe,skala}`,
`validation/si2-current-cpu-mpi-r2.log`, retained parity protocols in
`/tmp/cppaw-skala-mpi.Z4U801`. Restart SHA-256:
`22c48ba177e58ec0589cd7187e9c40540057c039075d4bd30c49467ee06baa7c`.
CPU serial executable for the applied warm probe:
`9e46ebcbb2d33eb2db04044a6e7dc851f867c7e215f04d463552841f23142ffd`.
CPU MPI executable:
`7425a9eb6c22fd35d335ef6599f46047752dc7c434e8bc8f00870fc51d820cd0`.
The CPU model hash is unchanged from the earlier checkpoint.

### Spark CPU parallel check and relaxation probes

Spark is used for CPU as well as GPU tests. Its updated
`nvhpc_gpu_fast_parallel` executable passes the same warm-start check with
one and eight MPI ranks, using `CPPAW_GPU_MODE=off`, `DEVICE='CPU'`, an empty
`CUDA_VISIBLE_DEVICES`, and one OpenMP/OpenBLAS thread per rank. The CPU
TorchScript model is the same export used on Terok. This is CPU execution of
a CUDA-capable NVHPC/NVPL binary, not a GPU-library-free executable; the latter
is tested separately by the GNU CPU-only build on Terok.

| Spark CPU 1/8-rank difference | Maximum absolute difference |
| --- | ---: |
| Total energy | 0 at the printed precision |
| Forces | 7.634e-10 H/bohr |
| Strain derivatives | 1.228e-9 H |
| Smooth-operator norms | 4.862e-8 |
| Electronic stationarity diagnostics | 2.290e-10 H |

Both runs have eight k points, 96 radial points, angular order 17, one
orientation and one shell. The eight-rank occupied maximum residual is
0.1405628149 H, so this is parallel consistency, not Skala convergence.
The mixed GNU/NVIDIA OpenMP warning remains visible and unresolved; a passing
numerical probe does not certify that runtime combination as safe. These
runs are not a CPU/GPU performance benchmark.

Evidence: Spark `validation/si2-current-spark-cpu-mpi8.log` and retained
protocols `/tmp/cppaw-skala-mpi.fH6p5H`. MPI executable SHA-256:
`9239704d5ddd10b718324476fea7fad50c10b49933bed2d0153f22869087fdb0`.
The PBE restart hash is the same as above; CPU model SHA-256:
`7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd`.

The same current executable also passes one/two-rank MPI parity with
`CPPAW_GPU_MODE=resident`, `DEVICE='CUDA'`, the CUDA model, and default
root-only Skala inference on Spark's single GPU. Both use the same PBE restart
and coarse grid as the CPU check:

| Spark GPU 1/2-rank difference | Maximum absolute difference |
| --- | ---: |
| Total energy | 0 at the printed precision |
| Forces | 7.981e-10 H/bohr |
| Strain derivatives | 1.980e-9 H |
| Smooth-operator norms | 1.235e-8 |
| Electronic stationarity diagnostics | 2.392e-10 H |

Evidence: `validation/si2-current-spark-gpu-mpi2.log` and
`/tmp/cppaw-skala-mpi.m4ECqE`. This closes the current-build GPU MPI parity
check, not multi-GPU/distributed-model validation, stationary-state finite
differences, or the mixed-OpenMP-runtime acceptance limitation. It ran alongside
a separate CPU crystal calculation and must not be used as a timing benchmark.

Comparing the one-rank CPU and GPU protocols above, with identical structure,
restart and quadrature, gives energy difference 1.539133e-7 H, force difference
1.157903e-9 H/bohr, strain-derivative difference 5.866493e-9 H,
smooth-operator norm difference 1.818030e-6, and electronic-diagnostic difference
3.572231e-9 H. These are observed backend differences for the float32 model,
not bitwise equivalence, stationary-state force accuracy, or a new reference
energy. The CPU and CUDA model exports have distinct hashes as recorded above.

Two bounded electronic-relaxation probes were also completed from that PBE
restart, with `DT=5`, `MPSI=100`, fixed nuclei/cell and applied Skala:

| Probe | Radial/angular grid | Steps | Final occupied RMS (H) | Final occupied maximum (H) | Final commutator (H) |
| --- | --- | ---: | ---: | ---: | ---: |
| Spark, CUDA model | 200 / 65 | 6 | 2.163144434e-2 | 4.118231448e-2 | 8.122257522e-3 |
| Terok, GNU CPU-only | 96 / 17 | 50 | 5.066137055e-3 | 2.203613351e-2 | 7.220628409e-3 |

The Spark residual was evaluated by a subsequent `CHECK=T`, `DT=1e-6` probe
of the final restart. It has overlap error 3.775e-15, Hamiltonian Hermiticity
error 4.164e-16 H, and integrated electrons 28.0000295156772 versus trace 28.
Its model XC energy is -42.1187507978530 H and total energy
-46.9375341367190 H. The maximum reported model-adjoint contraction error
is below 1e-14. No cell-stress finite difference was requested in this probe.

The CPU run records diagnostics throughout steps 300--349. Its RMS decreases
from 0.09426973264 to 0.005066137055 H; the final total energy is
-46.9384280625642 H. All 50 steps pass diagnostic-consistency checks, but the
final maximum residual and commutator remain far above 1e-6 H. Both probes
are **not converged**. Different grids and iteration counts prohibit treating
their energies as CPU/GPU parity data or quadrature-convergence evidence.

Evidence: Spark `validation/si2-relax-pilot-k8` and
`validation/si2-stationarity-gpu`; Terok `validation/si2-relax-coarse`.
The Spark diagnostic executable SHA-256 is
`787dae5e53cbfe0ba7a56070d8f19155d703767e47d3236d9921a79c652a653c`.
The pilot used the preceding Lebedev-order implementation; the subsequent
probe adds the read-only electronic diagnostics. The CPU executable/model
hashes are the same as in the preceding CPU MPI check.

### Positive-tau and setup-operator audit

The angular one-center **valence** tau now has an independent angularly exact
radial reference on the same domain. Their differences are of order 1e-14 H
for both Si atoms. Frozen-core tau is added separately in the production
reconstruction; it is not part of this diagnostic comparison.

On the warm PBE restart, atom 1 has a raw tau-minus-setup-kinetic difference
of 0.05111377703 H. Matching the outer radial domain contributes -0.03483613338 H
and the nominal AE ZORA factor contributes -0.006761831279 H, with a zero outer
surface term. The remaining gradient-form estimate minus the stored setup
matrix is 0.009515812379 H. No physical tau was changed to fit this number.

An explicit `!SPECIES!AUGMENT!GRID` sweep further separates the numerical
gradient and differential-operator discretizations from the stored HBS matrix:

| DMIN / DMAX (bohr) | Gradient minus numerical operator, atom 1 (H) | Numerical operator minus setup matrix (H) |
| --- | ---: | ---: |
| 1e-6 / 0.1 | 2.094133912e-5 | 9.494871040e-3 |
| 5e-7 / 0.05 | 5.219971007e-6 | 9.536840995e-3 |
| 2.5e-7 / 0.025 | 1.303286775e-6 | 9.465963259e-3 |

The first difference falls approximately quadratically with spacing. The second
does not disappear. These rebuilt setups also shift `R(NR-3)` and regenerate
partial waves; they are not a fixed-function quadrature test. HBS constructs
the pseudo waves inside a matching radius, retains its input tails outside,
and assigns `TPSPHI=(E-POT)*PSPHI`. Its stored kinetic action must therefore be
audited against that construction, not assumed equal to applying the nominal
nonrelativistic operator globally. This observation narrows the remaining
audit but does **not** certify the residual as harmless or as the cause of the
AlN EOS issue. Existing setup matrices and PBE reference energies are unchanged.

Evidence: Terok `validation/si2-tau-radial-half-r2` and
`validation/si2-tau-radial-quarter`. The first `si2-tau-radial-half` trial did
**not** refine the grid: historical `PARMS_STP` is not the current `AUGPARMS`
input, so the built-in setup was used. This trial is retained but excluded from
the sweep. The successful trials use explicit structure blocks and verify
the resolved grid in `si2.strc_out`.

### Geometry-cache checkpoint (2026-09-29)

The new bounded, per-rank partition cache stores only geometry derivatives.
It does not cache densities, neural-model adjoints, Hamiltonians or forces.
Exact atom/cell/grid, image-shell and force/stress-mode keys invalidate it;
insufficient budget or allocation failure selects the uncached computation.
The default is 256 MiB per MPI rank, with
`CPPAW_SKALA_PARTITION_CACHE_MB=0` retaining the uncached path. Array payloads,
not allocator overhead, count against the budget.

Model-free cache tests pass with local GNU, Spark NVHPC and Terok GNU builds.
They require bit-exact kernel reuse and test geometry/mode changes, parent
deallocation and exact, insufficient and zero budgets. The full local GNU
build/install also succeeds without Torch or NVIDIA linkage. The installer
now excludes standalone test executables and dangling optional module links.
The 19 Si2 and eight crystal Python tests pass.

Two-step applied-Skala comparisons start from identical PBE restarts, using
96 radial points, angular order 17, eight k points and a 20 Ry cutoff. They
compare both steps, not just the initial snapshot. Fixed-cell probes use
`DT=0.001`; stress probes use `DT=1e-6` and a very heavy moving cell.

| Cache-off/on comparison | Total energy (H) | Skala force (H/bohr) | Total force (H/bohr) | Largest smooth-operator norm difference |
| --- | ---: | ---: | ---: | ---: |
| Terok GNU CPU, fixed cell | 5.045e-13 | 3.660e-16 | 4.725e-9 | 1.400e-12 |
| Spark NVHPC resident GPU, fixed cell | 0 printed | 1.688e-9 | 1.489e-9 | 6.235e-8 |
| Terok GNU CPU, stress probe | 0 printed | 0 printed | 0 printed | 0 printed |
| Spark NVHPC resident GPU, stress probe | 0 printed | 7.351e-10 | 7.561e-10 | 2.406e-8 |
| Terok GNU CPU, eight MPI ranks, fixed cell | 2.203e-12 | 2.708e-15 | 2.699e-9 | 2.710e-11 |

The GPU stress probe differs by at most 1.677e-9 H in total strain derivatives
and 6.100e-13 H in the partition contribution. Both fixed-cell probes have
17,710 rows, zero first-step hits and 17,710 second-step hits, using 1,346,056
bytes. Stress probes allocate 4,463,016 bytes but correctly invalidate all
rows after tiny cell changes. They validate stress-path parity/invalidation,
not stress-path cache speedup.

The comparison uses a common CPU bound of 1e-10, an explicit GPU bound of
1e-7, and a separate total-force bound of 1e-8 H/bohr. An independent GNU
uncached repeat reproduces the same 4.725e-9 H/bohr total-force variation
after propagation, with Skala forces agreeing within 4e-16 H/bohr. Independent
GPU uncached repeats differ by up to 1.250e-8 in operator norms and 1.540e-9 H
in strain derivatives. Thus the end-to-end comparisons are not bit-exact
model evaluation, despite bit-exact reuse of the cached kernel values.
The original stricter failures remain in the evidence, not replaced by
fabricated passing results. The final comparator reanalysis passes the four
single-rank cases, and the eight-rank driver passes directly. Trials with
`DT=0.01` failed the overlap criterion and are excluded.

The two-step Spark fixed-cell wall times were 27.783 s uncached and 22.395 s
cached, including startup and the first cache fill. This is a single
integration-test pair, not a repeated steady-state performance benchmark.
Geometry caching primarily benefits fixed-geometry electronic iterations;
moving atoms/cells invalidate it.

The two-step test also exposed an unnecessary host update of `PROJ` in the
electronic stationarity diagnostic. Its device mapping could already have
been released, causing the second GPU step to abort. The diagnostic now uses
the host projection already produced by projection and its MPI reduction.
The wavefunction and Hamiltonian host updates remain unchanged.

Evidence in each host's `validation/`: Spark
`si2-cache-parity-fixed-gpu-dt0001`, `si2-cache-parity-gpu-stress-repeat`;
Terok `si2-cache-parity-fixed-cpu-dt0001` (including `off-repeat`) and
`si2-cache-parity-cpu-r2`, plus `si2-cache-parity-cpu-mpi8-r2` for eight ranks.
Their `provenance.json` files record inputs and
executables. The PBE restart, structure and CPU/CUDA models are unchanged.
Spark serial executable SHA-256:
`843efd00b00f3e88023445e1258e405e371272084f50fd4e13b0cf39f0408988`.
Terok GNU MPI executable SHA-256:
`7cb88d748357ad8f1d315b836b76c1a29f3cdc0537824cc6bf5b7edf433db43f`.
The eight-rank cache-off/on probe also has zero first-step hits, 17,710
second-step hits and 1,346,056 bytes summed over ranks, exercising distributed
atom ownership and the diagnostic reduction. This is not a CPU scaling
benchmark. The new cache has not yet been tested with multi-rank GPU execution;
the preceding GPU MPI parity evidence predates it.

### Cached relaxation and next performance targets

Spark completed another 100 applied-Skala steps, 349--448, from the preceding
50-step CPU coarse-grid restart, with the same 96/17 quadrature and `DT=5`.
Restart SHA-256: `c12b8273ab90b7dfccd764c9988db8dfee62e14778566f311417f883ce07a29f`.
Evidence: `validation/si2-relax-cached`, including initial-input hashes,
protocol and final restart. It used the Spark executable above, resident GPU
mode and one OpenMP/OpenBLAS thread. The first step filled 17,710 rows; the
remaining 99 steps each reused every row (1,753,290 hits, no further misses).

| Electronic diagnostic | Initial | Final |
| --- | ---: | ---: |
| Occupied RMS residual (H) | 5.066137077e-3 | 2.183459106e-5 |
| Occupied maximum residual (H) | 2.203613406e-2 | 3.445278785e-5 |
| Occupation commutator (H) | 7.220628871e-3 | 3.399011398e-5 |

All 100 steps pass diagnostic consistency, with maximum overlap error
1.228e-9 and Hamiltonian Hermiticity error 1.640e-15 H. Final total energy is
-46.9388733875965 H. The final five-step window **fails** the 1e-6 H residual
and commutator acceptance target; no stationary-state acceptance is claimed.
The mixed-OpenMP-runtime warning remains unsuppressed. NVHPC also reports
IEEE invalid/divide-by-zero/underflow/inexact flags at exit despite finite
reported diagnostics; their origin has not been isolated or certified harmless.

The run-time report gives 555.6 s elapsed and the following accumulated
**CPU-time** breakdown: reconstruction 187.2 s (34%), adjoints 351.6 s (63%),
model call 3.8 s (1%). These are not CUDA-event model timings. The separate
accelerator wall-clock profile also puts the atom-block phase first:
`PAW_ETOT_SPHERE` 557.357 s versus `PHASE_ETOT_WAVES` 566.946 s. Nested timers
must not be added. This identifies reconstruction/adjoint work as the current
priority, not another FFT/library substitution for this small Si2 case.

Source inspection identifies the next opportunities, not yet implemented:

1. With CUDA and default `DISTRIBUTEGPU=F`, `WAVES$DENMAT` assigns all atom
   blocks, including CPU reconstruction and adjoints, to MPI rank 1. CPU mode
   distributes whole atoms, so Si2 exposes at most two atom-block workers.
   Point-level parallelism can use more cores without duplicating the GPU model.
2. Hoist cell inversion and image-bound calculations; cache interpolation and
   orbital basis/gradient data by geometry with an explicit memory budget.
   Dense orbital-pair kernel caches would scale poorly and are not the target.
3. Batch reconstruction and adjoint contractions on the GPU, retaining basis
   data and iteration-local density matrices and accumulators there. Shared
   Hamiltonian, force and stress writes need explicit, tested reductions.
4. Avoid unused Hessians/Jacobians in forward-only paths. Do not remove needed
   force/stress terms, change quadrature to improve timings, or split a complete
   coupled atom-block model into independent pointwise model evaluations.

These optimizations do not close the physical acceptance gates below.

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
