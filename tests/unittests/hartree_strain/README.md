# Hartree strain derivatives

This test links the production `POTENTIAL_HARTREE` kernel, reciprocal-space
structure factors and solid spherical harmonics. It compares all six symmetric
strain derivatives against centered energy differences at fixed fractional
nuclei, fixed multipole moments and fixed volume-scaled density coefficients.
Analytic Gaussian radial factors include their reciprocal-radius derivatives
and inverse-volume normalization.

Three cases use equal multipole limits, mixed monopole/quadrupole limits and
mixed dipole/quadrupole limits. Each also permutes three atoms and checks that
energy and the full strain tensor do not depend on atom order. This detects
using the previous atom's species when constructing angular strain terms.

Run against an existing serial build:

```sh
make -C bin/Build_fast -f Makefile \
  -f ../../tests/unittests/hartree_strain/driver.mk hartree-strain-test
```

The test has no Skala or CUDA requirement. It verifies the isolated discrete
electrostatic kernel, not converged material stresses or a full PAW calculation.
