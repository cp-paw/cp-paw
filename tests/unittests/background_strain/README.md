# One-center background strain test

Run against a configured build with

```sh
make -C bin/Build_fast -f Makefile \
  -f ../../tests/unittests/background_strain/driver.mk background-strain-test
```

The test calls the production all-electron and pseudo one-center Hartree
functions. Their separately returned background energies must equal the
change from a zero-background evaluation. They remain included in the total
Hartree energy and must not be added a second time.

At fixed compensating charge, its density varies as inverse cell volume.
The explicit strain derivative of the one-center energy difference is
`-(E_background_AE - E_background_PS)` on each diagonal, with zero shear.
Finite differences cover a single diagonal deformation, uniform dilation,
and symmetric shear for zero, positive, negative, and small background
densities. The integration grid, local densities, and multipoles stay fixed.
This kernel test is not a complete stationary-state stress validation.
