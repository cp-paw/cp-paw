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
