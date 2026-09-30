# Joint PAW reconstruction and periodic quadrature

Run the compiler-independent source, grid and partition tests against an
existing CP-PAW build:

```sh
bash tests/unittests/skala_reconstruction/run.sh bin/Build_skala_cpu_fast
```

The periodic partition tests cover atom and point derivatives, all nine
affine cell derivatives with fixed Cartesian local quadrature offsets,
self-image descriptors, skew cells, integer lattice translations, and
changes of the nearest-image enumeration. The sum over all atom/image
weights must be one. These tests do not evaluate the neural model.

The same suite tests the Lebedev rules through order 65, including analytic
Cartesian moments and rounding a requested exactness of 64 to the 1454-point
rule. See [LEBEDEV.md](LEBEDEV.md) for coefficient provenance and test coverage.

The partition-cache tests compare stored derivatives bit for bit with fresh
kernel evaluations, including skew cells and self-image strain derivatives.
They cover changes to atom positions, cell, grid coordinates and size, source
atom count, target atom, image shell and force/stress mode, and verify parent
deallocation and exact/insufficient/zero storage budgets. Invalid cache-limit
environment values must fail explicitly. The cache stores only geometry; it
does not validate the neural model or electronic convergence.

The source-geometry cache has independent tests for AE/pseudo partial-wave
values, gradients, Hessians and frozen-core fields at periodic source images.
These compare cached and direct fields, density-matrix adjoints, forces and
strain moments in both collinear spin modes. Changed density matrices and model
adjoints must not invalidate geometry or reuse electronic outputs. Changes in
geometry, row count, setup arrays and derivative mode must invalidate it.
Tests include empty source support, full and partial row coverage, exact and
insufficient budgets, disabled caching and invalid environment settings.

The source contractions first form density-matrix products with the partial
waves and their three gradients. This moves the Hessian contraction out of
the partial-wave pair loop without changing the forward fields or reverse
map. An independent direct-pair reference checks 1, 6 and 11 partial waves,
both collinear spin modes, nonsymmetric complex density matrices, and
accumulation into existing complex adjoints. It compares fields, coordinate
derivatives and density-matrix contractions to an absolute bound of 1e-12,
including calls that do not request Hessians. The complete reconstruction
suite also runs in the CPU-only Skala CI job.

## Smooth native-grid interpolation

The native-grid interpolant blends two nine-node, degree-eight Lagrange
stencils centered on the endpoints of each grid interval. For interval
coordinate `t`, the blend is `s(t)=t^3 (10-15t+6t^2)` and the interpolant is
`(1-s) P_left + s P_right`. Both polynomials reproduce degree eight.
The blend and its first two endpoint derivatives select the same centered
polynomial on either side of a grid node. The piecewise interpolant is C2.
The union has ten nodes per direction with periodic wrapping.

Point derivatives differentiate both the polynomials and the blend. The
reverse map scatters with the identical tensor-product weights. CPU and
OpenACC backprojection call the same weight routine. Tests cover polynomial
reproduction, gradient finite differences and second-derivative continuity
across internal and periodic faces. The grid adjoint tests include points
on stencil boundaries as well as generic points and skew cells.

The earlier single eight-node Lagrange stencil was only C0 across interval
changes. This is a revised interpolation discretization, not a bitwise
equivalent optimization. Electronic states and finite differences must be
reconverged when comparing the two. Smoothness does not imply exact
translation covariance, density-grid convergence or physical force accuracy.

## Periodized Becke weights

The previous common finite cluster gave its edge images inequivalent
neighbour environments. Normalization inside each cluster alone did not
provide a periodic partition of unity for the collection of atom grids.
Increasing `IMAGESHELLS` reduced this error but did not remove it at the
default shell count.

The new kernel constructs a translated seed for each atom image. With
`f = inverse(H) x`, the compact window is a product of three factors that
are one for `abs(f[d]) <= S`, zero for `abs(f[d]) >= S+1`, and the C2
quintic `1 - 10 z^3 + 15 z^4 - 6 z^5` in between. `S` is `IMAGESHELLS`.
The ordinary triple-cubic Becke pair switch is `s_AB(x)`. Define

```text
q_A(x) = C(x) product_Bn [1 - C(R_Bn-R_A) (1-s_ABn(x))]
w_A0(r) = q_A(r-R_A) / sum_Bn q_B(r-R_Bn).
```

The denominator includes every nonzero compact seed. Each seed has its own
translated neighbour environment. Both seed support and neighbour membership
are tapered, so enumeration changes add or remove only zero contributions.
The reverse includes the explicit cell dependence of both windows, not just
the motion of image centres. No density rescaling or fitted quadrature-volume
correction is used.

This is a revised finite-shell discretization, not an algebraically identical
implementation of the old finite-cluster rule. Shell convergence of the
energy and the self-image descriptor remains necessary even though the
periodic weights already sum to one. Source reconstruction has its separate
radial-support image search and must not be truncated to these neighbours.

## Model-free integration probe

```sh
make -C bin/Build_skala_cpu_fast -f Makefile \
  -f ../../tests/unittests/skala_reconstruction/driver.mk skala-partition-probe
bin/Build_skala_cpu_fast/unit-tests/partition_measure.x 1 200 53 1 1e-4
```

Arguments are image shells, radial points, Lebedev exactness, orientations,
and an optional tolerance. The Si2 primitive cell has volume 270.011394
bohr^3. The probe integrates a constant and the first reciprocal cosine mode;
the latter must vanish. An explicit tolerance tests both errors relative to
the cell volume. Without a tolerance, successful execution is not a
convergence certificate.

Constant integration, electron-number convergence and model-energy
convergence are distinct tests. In particular, plane-wave cutoff and native
grid interpolation cannot affect this model-free geometry probe.

## GPU Source Reverse

The source reverse test compares a 257-row batch with the independent row
path, crossing the internal tile boundary. It covers one/two spin channels,
complex nonsymmetric input and pre-existing matrices, force omission,
source forces and image moments, and full/partial/disabled host caches.
CPU-only builds exercise the fallback. To require actual device coverage
with an OpenACC build on an NVIDIA GPU, run the compiled test with

```sh
CPPAW_SKALA_SOURCE_BACK_ACC=1 CPPAW_SKALA_SOURCE_BACK_ACC_MIN_ROWS=1 \
CPPAW_SKALA_TEST_REQUIRE_SOURCE_ACC=1 ./unit-tests/skala_reconstruction.x
```

Also run without the coverage requirement and with
`CPPAW_SKALA_SOURCE_BACK_ACC_MB=0` to exercise device-budget fallback.
