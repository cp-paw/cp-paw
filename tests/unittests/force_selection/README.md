# Nuclear-force selection

Run `python3 tests/unittests/force_selection/input_checks.py /path/to/paw.x`.
Use `--mpi-ranks 2 --mpiexec /path/to/mpirun` for an MPI build and
`--output /path/to/new-directory` to retain raw results.

The test creates a shared conventional Si2 restart, then compares omitted,
explicitly enabled and explicitly disabled `!PSIDYN FORCE` inputs.
Electronic wave-function and Lambda restart payloads must agree within
`1e-10`, with identical metadata. The numeric cell stored in the electronic
restart may differ by at most `1e-12` Bohr in moving-cell comparisons, since
even the heavy reference cell can develop nonzero entries of order `1e-28`.
All fixed-cell comparisons require an identical stored cell. Energy trajectories must agree within
`1e-10` Ha. The numerical comparison permits floating-point variation between
independently initialized FFT/linear-algebra backends. The fixed-geometry
`FORCE=F` calculation must not write force-trajectory records.

Matched atom- and cell-dynamics controls verify that `FORCE=F` cannot disable
required forces, including explicit `STRESS=F` during cell motion. Fixed-cell
`STRESS=T` controls likewise require complete nuclear forces and unchanged
electronic evolution. Their force trajectories must agree within `1e-10` Ha/bohr
with `FORCE=T`, with identical step/time metadata. This is an implementation
regression, not a physical force or
stress convergence test. QM/MM also forces nuclear-force work on in the
production selector, but is not exercised by this Si2 fixture.
An additional displaced-atom control requires a nonzero force trajectory
and verifies that its position record is the evaluated initial geometry,
before propagation. This guards against sampling after `ATOMS$SWITCH`, which
advances positions and clears forces. The final-step omission is unchanged.

`cell_residency.py /path/to/profile/paw.x --output /path/to/new-directory`
requires a GPU-enabled profile build. It compares four moving-cell steps
with GPU mode off, transfer and resident from one conventional Si2 restart.
The cell must move measurably. Cell histories, wave/Lambda payloads, energies
and forces must agree within `1e-8` (electronic cell metadata within `1e-9`).
The second atom is displaced by `0.01` Bohr along x in each restart and atom
dynamics is enabled. The reference must have nonzero trajectory forces,
not symmetry-enforced zeros or zeros caused by sampling after force reset.
The atomic position histories must also agree within `1e-8` Bohr.
Profile counters must prove that setup left PSIM resident and host propagation
explicitly retrieved it. This catches a stale device copy masking the host
update when cell motion selects CPU propagation. Use `--mpi-ranks` for MPI.
