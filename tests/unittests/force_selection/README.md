# Nuclear-force selection

Run `python3 tests/unittests/force_selection/input_checks.py /path/to/paw.x`.
Use `--mpi-ranks 2 --mpiexec /path/to/mpirun` for an MPI build and
`--output /path/to/new-directory` to retain raw results.

The test creates a shared conventional Si2 restart, then compares omitted,
explicitly enabled and explicitly disabled `!PSIDYN FORCE` inputs.
Electronic wave-function and Lambda restart payloads must agree within
`1e-10`, with identical metadata. Energy trajectories must agree within
`1e-10` Ha. The numerical comparison permits floating-point variation between
independently initialized FFT/linear-algebra backends. The fixed-geometry
`FORCE=F` calculation must not write force-trajectory records.

Matched atom- and cell-dynamics controls verify that `FORCE=F` cannot disable
required forces. Their force trajectories must agree within `1e-10` Ha/bohr
with `FORCE=T`, with identical step/time metadata. This is an implementation
regression, not a physical force or
stress convergence test. QM/MM also forces nuclear-force work on in the
production selector, but is not exercised by this Si2 fixture.
