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

## Partition Derivative Cache

`cache_parity.py` runs two applied-functional steps from the same restart with
the geometry cache off and on. It requires exact row reuse on the second
cached step and compares both steps' total/XC energies, total and partition
forces, operator norms and electronic residuals at fixed nuclei and cell.
The optional `--stress` probe also evaluates stress using a very heavy moving
cell. Even tiny cell changes correctly invalidate the cache, so this variant
checks parity and hit/miss accounting without requiring second-step hits.
The default 96/17 grid is only an integration test, not converged quadrature.

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
