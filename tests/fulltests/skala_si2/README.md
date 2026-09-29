# Periodic Si2 checks

`skala_si2.strc` is the two-atom primitive diamond cell, not an isolated Si2
molecule. It uses eight k points. `run_isolated.sh` performs a cold one-step
integration test in a fresh directory; this is not an SCF convergence test.

`RADIALPOINTS` selects the model atom-grid radial rule; `LEBEDEVEXACTNESS`
selects the minimum angular polynomial order. The maximum is 65, with 1454
directions; a request for 64 selects 65. The default remains 53. Orientations
average rotated rules without increasing their polynomial exactness.

## Electronic stationarity

Enable `!DFT!SKALA CHECK=T` and inspect the completed protocol:

```sh
python3 tests/fulltests/skala_si2/stationarity.py /path/to/si2.prot
python3 tests/fulltests/skala_si2/stationarity.py /path/to/si2.prot \
  --residual-tolerance 1e-6 --commutator-tolerance 1e-6 --last 3
python3 -m unittest discover -s tests/fulltests/skala_si2 -p test_stationarity.py
```

Without the two explicit convergence tolerances, the script checks diagnostic
consistency only. The convergence mode limits the **maximum occupied residual**
and occupation commutator, in Hartree, over the requested final window; the RMS
is also reported. Overlap and Hermiticity limits apply to every recorded step.
The example tolerances are not universal force/stress accuracy guarantees.
Neither mode certifies the global electronic minimum or grid convergence.
`APPLY=F` measures the conventional-XC Hamiltonian; `APPLY=T` measures Skala.

`CPPAW_SKALA_SCF_DETAIL=1` additionally reports each global k-point, spin and
band, its weighted occupation, residual norm, Hamiltonian expectation and
maximum occupation-commutator element. Use `stationarity.py --bands` to read
these records and verify their weighted RMS and maxima against the aggregate
diagnostics. Empty-state residuals are reported but do not enter the occupied
RMS or maximum. This opt-in diagnostic does not change propagation.

A PBE restart must have the same cell, positions, k mesh, band layout and
Fourier basis as the Skala probe. A small time step preserves that starting
state but cannot make it stationary for a different functional. A one-k-point
restart cannot be used for this eight-k-point test.

For fixed-occupation relaxation, test `!PSIDYN SAFEORTHO=T` explicitly as
well as the current applied-Skala default `F`. The former uses the
occupation-weighted constraint update and does not require diagonalizing
the Hamiltonian within an equally occupied subspace. This distinction does
not remove the occupied--empty commutator requirement. Dynamical occupations
through `!MERMIN` still require `SAFEORTHO=F`.

When comparing fictitious masses, hold `MPSICG2` fixed to isolate the overall
mass from the reciprocal-space preconditioner. For example, at `DT=5`,
the automatic coefficient for `MPSI=100` is `0.3166286988823056`. Specify that
coefficient explicitly when testing `MPSI=25`. Use the same restart, geometry,
grid, friction and stationary window, and retain failed checks. These are
controlled relaxation settings, not universal cold-start recommendations.

## PAW setup grid

The partial-wave setup grid is independent of the Skala atom-grid quadrature.
The built-in `SI_.75_6.0` uses `DMIN=1.E-6 DMAX=0.1 RMAX=20.`. To refine it,
provide a complete `!SPECIES!AUGMENT` block, including its `!GRID` block, or use
the current external `AUGPARMS` format (`!ACNTL!AUGMENT`). Check the resolved
`*.strc_out` file to verify which settings were actually used.

The historical `PARMS_STP` file attachment in older test controls does **not**
override the built-in setup. Changing that `stp.cntl` alone therefore does not
constitute a radial convergence test. Change the resolved `AUGMENT!GRID`
parameters for a setup-grid convergence test.

## Orbital Energy Derivative

`orbital_fd.py` checks a fixed-geometry, fixed-occupation direction in orbital
space. Every single-step calculation starts from the same restart, applies
the initial PAW orthonormalization, and then rotates two bands before rebuilding
projections and densities. A nonzero center angle provides a measurable signal
even when the original restart is nearly stationary. The weighted occupations
already include the k-point weight and spin multiplicity.

```sh
python3 tests/fulltests/skala_si2/orbital_fd.py \
  --executable /path/to/paw.x --model /path/to/model.fun \
  --restart /path/to/si2.rstrt --structure /path/to/si2.strc \
  --output /path/to/new-orbital-test-directory \
  --kpoint 1 --bands 4 5 --center-angle 0.02 --steps 0.01 0.003 0.001
```

The driver repeats the center calculation and records central energy
differences and the Hamiltonian derivative separately for each step size.
Without `--absolute-tolerance`, it reports measurements, not a passed test.
With that option every step must meet the specified absolute bound plus
`--relative-tolerance` times the analytic derivative magnitude. Establish
truncation and model-precision errors before interpreting those bounds.
The default 96/17 quadrature is not a physical convergence claim.

Use `--reference-pbe` to retain the conventional PBE Hamiltonian (`APPLY=F`).
For GPU checks select `--device CUDA --gpu-mode resident` with a compatible
model. Complex, unpacked orbitals also allow `--mode IMAG`. Imaginary rotations
are rejected for the packed real representation. The underlying diagnostic
is opt-in through `CPPAW_SKALA_ORBITAL_ROTATION='K SPIN I J ANGLE REAL|IMAG'`
and requires Skala `CHECK=T`. It is intended for a single energy evaluation,
not production propagation. This test neither converges the electrons nor
validates stationary-state ionic forces or stress.

To exercise complex Bloch orbitals, create a cold PBE seed with the reduced
`3 1 1` mesh. Its second irreducible k point is `(1/3,0,0)` and is unpacked.
It deliberately contains high-energy components and is not an SCF reference.

```sh
python3 tests/fulltests/skala_si2/orbital_seed.py \
  --executable /path/to/paw.x --output /path/to/new-seed
python3 tests/fulltests/skala_si2/orbital_fd.py \
  --executable /path/to/paw.x --model /path/to/model.fun \
  --restart /path/to/new-seed/si2.rstrt --structure /path/to/new-seed/si2.strc \
  --output /path/to/new-complex-real-test --kpoint 2 --mode REAL \
  --steps 0.01 0.003 --absolute-tolerance 1e-5 --relative-tolerance 1e-3
```

Repeat with `--mode IMAG` and a fresh output directory, and use
`--reference-pbe` for the host control. The coarse-grid tolerances above
test integrated energy/operator consistency, not physical accuracy. This
probe is sensitive to a missing density-space Fourier projector in the
scalar adjoint. The forward density is band-limited before Skala evaluation,
so its reverse scalar potential must use the same projector. Positive tau
is formed directly on the native grid and must not inherit that filter.

Use `--direction kinetic --band 4` instead of `--bands 4 5` to test outside
the current band space. From the common orthonormalized restart the diagnostic
forms `chi = |G+k|^2 psi_i`, removes its PAW-metric projection onto every
current band twice, and normalizes it with the full PAW overlap operator.
It then evaluates `psi_i(theta) = cos(theta) psi_i + sin(theta) chi` and the
analytic derivative `2 f_i Re <d psi_i/d theta | H psi_i(theta)>`. This is
not another occupied--empty rotation within the computed band space.
`--mode IMAG` starts from `i |G+k|^2 psi_i` and requires unpacked orbitals.
The direction norm, unit-metric error, and orthogonality error are reported
and checked, including a common direction norm across the finite difference.

The corresponding runtime diagnostic is
`CPPAW_SKALA_ORBITAL_TANGENT='K SPIN BAND ANGLE REAL|IMAG'`. It requires
fixed occupations, collinear orbitals, and `CHECK=T`, is mutually exclusive
with `CPPAW_SKALA_ORBITAL_ROTATION`, and reports only the first energy
evaluation. The driver still requires exactly one completed step per leg.
Use a PBE reference and more than one finite-difference step before interpreting
the Skala result. Passing one external direction does not establish the
entire gradient, electronic stationarity, or a stationary force derivative.

## Model Precision Diagnosis

`model_precision.py` creates an isolated Float64 **diagnostic copy** of the
hash-pinned Skala-1.1-rev1 CPU or CUDA TorchScript export. It does not replace
the production model, change CP-PAW defaults, retrain weights, or supply an
official higher-precision model. The published Float32 weights are promoted
exactly. Explicit Float32 casts and factories in the scripted methods are
rewritten only at their dtype operands, preserving shared shape/index values.
Protocol metadata is preserved and an additional diagnostic marker is added.

```sh
python3 tests/fulltests/skala_si2/model_precision.py \
  /path/to/skala-1.1-rev1.fun /path/to/new-precision-check.fun --device cpu
python3 -m unittest discover -s tests/fulltests/skala_si2 -p test_model_precision.py
```

For CUDA, use the pinned CUDA export and `--device cuda`. PyTorch with the
selected backend is required. The graph inspection uses private TorchScript
APIs and is deliberately limited to these exports, not arbitrary future
models. It checks parameter/buffer values, dtype inventories, serialization,
metadata and synthetic energy derivatives for seven feature groups. Both
fine finite-difference steps (1e-5 and 1e-6) must have absolute errors at most
1e-8. A coarse 1e-4 step above that bound must improve at least twentyfold
at 1e-5. This allows truncation error without accepting nonfinite results.

Only a validated model and its hash-bearing JSON report are published, and
existing output files are never overwritten. Use that `.fun` file explicitly
with `orbital_fd.py --model` to separate arithmetic error from the integrated
PAW derivative. Sampled derivative agreement does not establish electronic
stationarity, force/stress accuracy, or an appropriate production precision.

## Geometry Caches

`cache_parity.py` runs two applied-functional steps from the same restart with
the geometry cache off and on. It requires exact row reuse on the second
cached step and compares both steps' total/XC energies, total and partition
forces, operator norms and electronic residuals at fixed nuclei and cell.
The optional `--stress` probe also evaluates stress using a very heavy moving
cell. Even tiny cell changes correctly invalidate the cache, so this variant
checks parity and hit/miss accounting without requiring second-step hits.
The default 96/17 grid is only an integration test, not converged quadrature.
`--cache source` instead checks the source partial-wave/core geometry cache.
The other cache is disabled in both legs to isolate the selected cache.
Source entries are prepared before the forward pass, so hits are expected
already on the first step, including forward and reverse evaluations of both
Si source atoms. The budget can cover only part of a block. This probe requires
the same nonzero coverage on both fixed-geometry steps and checks all remaining
rows are counted as misses. The model-free reconstruction test separately verifies exact
geometry reuse, changed density matrices, setup/geometry/mode invalidation
and memory-budget fallback. Cached source geometry currently stays on the host.
`--cache-mib` changes the enabled-leg budget from its 256 MiB default.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CPPAW_GPU_MODE=resident \
python3 tests/fulltests/skala_si2/cache_parity.py \
  --executable bin/nvhpc_gpu_profile/paw_nvhpc_gpu_profile.x \
  --model /path/to/skala-1.1-rev1-cuda.fun --device CUDA --tolerance 1e-7 \
  --restart /path/to/si2.rstrt --structure /path/to/si2.strc \
  --output /path/to/new-validation-directory
```

For CPU-only validation, use a GNU Skala executable, CPU model and
`--device CPU`; `--mpi-ranks N --mpiexec /path/to/mpirun` is optional.
The output directory must not exist. Inputs and executable hashes, environment,
protocols, comparison results and wall times are retained. The timings include
startup and the first uncached step and are not steady-state benchmark results.
`CHECK=T` now reports atom-block and operator diagnostics on every step.
The default timestep is 0.001: unlike a one-step snapshot, a two-step test also
exercises propagation and the subsequent constraint forces, for which an
extremely small timestep can amplify roundoff. `--timestep` overrides it.
The default comparison tolerance is 1e-10. The CUDA example allows 1e-7 for
float32 model-adjoint variability. This is distinct from the bit-exact cache
kernel test and does not establish physical force or stress accuracy.
Total forces have a separate default bound of 1e-8 H/bohr, adjustable with
`--total-force-tolerance`. The stricter common bound still applies to the
Skala and partition force contributions and all other diagnostics. Independent
uncached repeats can distinguish propagation roundoff from cache effects.
