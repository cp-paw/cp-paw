# Reciprocal Setup Interpolation

This model-free test calls `SETUP$GETFOFG` for all five supported radial
quantities using an exactly interpolable cubic. It checks values and radial
derivatives against analytic expressions, scalar calls, and a reversed vector
order. The inputs include identical radii, a chain of distinct radii separated
by less than the former absolute reuse threshold, and widely separated radii.

Finite differences check isotropic reciprocal scaling with its volume factor
and a perturbation of one radius independently of its neighbours. Approximate
reuse of a previous value fails the latter derivative check.

Run against an existing CPU or NVIDIA build with

```sh
make -C bin/Build_fast -f Makefile \
  -f ../../tests/unittests/setup_reciprocal/driver.mk setup-reciprocal-test
```

This tests the interpolation interface, not the accuracy of the setup transform
or full-system force and stress convergence.
