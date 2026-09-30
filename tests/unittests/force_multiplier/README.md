# Occupation-weighted PAW force multiplier

Run from a configured serial build with either GNU Fortran or NVHPC.

```sh
make -C bin/Build_fast -f Makefile \
  -f ../../tests/unittests/force_multiplier/driver.mk force-multiplier-test
```

The equation-of-motion matrix `RLAM0` is not generally Hermitian.
The constraint multiplier in the real PAW Lagrangian is the Hermitian part
of `RLAM0 * diag(OCC)`. The force contraction therefore uses minus that
matrix, not `-RLAM0(i,j) * (OCC(i)+OCC(j))/2`.

The test differentiates an independently constructed projector Lagrangian
containing both the Hamiltonian and overlap-augmentation terms. It compares
four complex directional derivatives against central finite differences for
equal, fractional, partially empty and entirely empty occupations. A fifth
case tests a non-Hermitian perturbation of the weighted multiplier, for which
only its Hermitian part contributes to the real Lagrangian. Equal occupations
retain the previous result. Unequal occupations explicitly expose the old
formula's error, including arbitrary equation-of-motion entries in an empty
column. This algebraic check complements integrated force and strain tests.
The actual `WAVES_ADDOPROJ` contraction is checked against the same gradient
in complex storage and against an independent real-matrix result in packed
time-reversal storage.
