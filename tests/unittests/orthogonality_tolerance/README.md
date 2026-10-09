# PAW constraint tolerance

After building CP-PAW, run

```sh
make -C bin/Build_fast -f Makefile \
  -f ../../tests/unittests/orthogonality_tolerance/driver.mk orthogonality-tolerance-test
python3 tests/unittests/orthogonality_tolerance/input_checks.py bin/Build_fast/paw.x
```

This links the production real and complex occupation-weighted constraint
solvers. The test checks the unchanged default of 1e-8, an overlap defect
that is accepted at that default but corrected at 1e-12, and independent
overlaps after solving at 1e-8, 1e-12 and 1e-14. The matrix cases include
unequal and zero occupations. Nonfinite and out-of-range setters must fail.
The overlap stopping criterion is distinct from electronic stationarity
and from the accuracy of physical forces.
The parser checks run isolated conventional Si2 inputs, including both
orthogonalization modes, the default, an explicit tight tolerance, and
rejection of incompatible or out-of-range settings.
With an MPI executable, add `--mpi-ranks 2` to the Python command. Both
serial and two-rank parser checks run in the conventional CPU CI job.

For an NVIDIA GPU build, the same executable can also be run with
`CPPAW_GPU_MODE=resident CPPAW_CUBLAS_ACC_ORTHO_X_RESIDENCY=1` to exercise
the resident real solver. The complex solver uses its existing path.
