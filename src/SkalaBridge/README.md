# Skala 1.1 FTorch bridge

This optional shared library isolates LibTorch and C++ from CP-PAW's legacy
Fortran build. It uses FTorch v1.0.0 for array-to-tensor transfers and a small
LibTorch protocol adapter for Skala's dictionary input and autograd output.
The setup path builds the pinned FTorch source even if another installation is
visible. CMake users can opt into an external package with
`-DCPPAW_SKALA_USE_SYSTEM_FTORCH=ON`.

The scientific input contract is deliberately PAW-specific. Host arrays use
`density(npoint,2)`, `grad(npoint,3,2)`, and `kin(npoint,2)`. The bridge converts
them to Skala protocol-v2 tensors and returns derivatives in the host layout.
`kin` must be the positive kinetic-energy density
`0.5*sum_i |grad psi_i|^2`; CP-PAW's existing Laplacian-gauge `RHOKIN` is not a
drop-in replacement.

The native PAW call path is being integrated in stages. At this point the
`!CONTROL!DFT!SKALA` block assembles complete hybrid PAW atom blocks, evaluates
the model, and maps its density, density-gradient, and positive-tau adjoints
back to the one-center density matrices and smooth wave-function grid. The
smooth scalar operator includes the Fourier-space divergence of the gradient
adjoint. These operators are validated diagnostically; they do not replace
CP-PAW's conventional XC energy or Hamiltonian by default. `APPLY=T` enables
the experimental electronic operator, including the generalized Kohn-Sham
positive-tau term. The smooth scalar, positive-tau, and one-center parts can be
isolated with their component switches. Its safe default is `F`. Experimental
analytic atomic forces include the explicit model-coordinate, moving local
grid, interpolated-primitive, and smooth/local partition terms in addition to
the PAW projector response. The smooth scalar adjoint is also contracted with
the translated pseudo-core densities after Skala inference; this contribution
cannot be taken from the earlier conventional-XC potential. Analytic stress is
available when CP-PAW requests cell stress. It includes the model-coordinate,
quadrature-volume, hybrid-partition, radial-blend, interpolated smooth-field,
pseudo-core, positive-tau, and existing PAW projector responses.

```text
!SKALA
 MODEL='path/to/skala-1.1-rev1-cuda.fun'
 DEVICE='AUTO'
 RADIALPOINTS=100
 LEBEDEVEXACTNESS=17
 LEBEDEVORIENTATIONS=1
 APPLY=F
 APPLYSMOOTH=T
 APPLYTAU=T
 APPLYONECENTER=T
 INTERPOLATEPARTITION=F
 CHECK=F
!END
```

`CHECK=T` performs a one-time central finite-difference check of the model
adjoint and its PAW one-center density-matrix contraction. It also verifies
the discrete integration-by-parts identity for the density-gradient adjoint
and compares the occupied-state expectation of the positive-tau Hamiltonian
with the primitive `integral v_tau*tau`. `CHECK=F` omits the diagnostic field
copies and the additional wave-function overlap. Force diagnostics report the
total energy at full double precision and list the pseudo-core contribution
separately, so central differences can be evaluated without the rounded energy
summary. Stress diagnostics report each Skala contribution and the complete
`D E / D STRAIN` tensor before CP-PAW converts it to its internal stress sign
convention.

`APPLY=F` leaves CP-PAW's conventional XC functional active. Consequently,
subtracting forces from otherwise identical `APPLY=T` and `APPLY=F` runs gives
the Skala-minus-conventional-XC force, not the derivative of the Skala model
energy alone. End-to-end force finite differences must compare total energies
and analytic forces from the same applied functional at an electronically
converged state.

End-to-end stress differences likewise require an electronically converged
restart. Apply `F = I + strain` to the restart cell and atom positions while
retaining the old reciprocal basis stored with the wave functions, then compare
`(E(+h)-E(-h))/(2h)` with the reported `D E / D STRAIN`. Local PAW
radial/Lebedev offsets remain fixed in Cartesian space under this deformation;
only their atom centers follow the affine cell motion. The interpolated smooth
fields therefore contribute through the negative local-offset response rather
than an affine deformation of the radial grid. The Si2 validation covers
isotropic, uniaxial, and symmetric-shear strains and 1/4-rank MPI parity.

Given an electronically converged restart and its matching structure, the
force and stress checks can be repeated with:

```sh
cd tests/fulltests/skala_si2
PAWX=/path/to/ppaw SKALA_MODEL=/path/to/model.fun \
SKALA_RESTART=/path/to/si2.rstrt SKALA_STRUCTURE=/path/to/si2.strc \
  ./force_fd.sh 2 1

PAWX=/path/to/ppaw SKALA_MODEL=/path/to/model.fun \
SKALA_RESTART=/path/to/si2.rstrt SKALA_STRUCTURE=/path/to/si2.strc \
  ./stress_fd.sh isotropic

PAWX=/path/to/ppaw SKALA_MODEL=/path/to/model.fun \
SKALA_RESTART=/path/to/si2.rstrt SKALA_STRUCTURE=/path/to/si2.strc \
MPI_RANKS=4 ./mpi_parity.sh

PAWX=/path/to/ppaw SKALA_MODEL=/path/to/model.fun \
SKALA_RESTART=/path/to/si2.rstrt SKALA_STRUCTURE=/path/to/si2.strc \
  ./quadrature_convergence.sh
```

The force arguments select the one-based atom and Cartesian axis. The other
stress cases are `xx` and `xy`. `MPI_RANKS` selects an MPI run; both drivers
default to a step of `3e-4` and an absolute tolerance of `2e-3`.
`SKALA_RADIAL_POINTS`, `SKALA_LEBEDEV_EXACTNESS`, and
`SKALA_LEBEDEV_ORIENTATIONS` override the local-grid settings in the force,
stress, and MPI drivers for quadrature-convergence checks.
`quadrature_convergence.sh` defaults to radial counts `100 200 400` at Lebedev
exactness 17 and reports every result relative to the final, finest case. The
lists can be changed with `SKALA_RADIAL_POINT_LIST` and
`SKALA_LEBEDEV_EXACTNESS_LIST`. `SKALA_LEBEDEV_ORIENTATION_LIST` adds a
convergence sweep over deterministic, equally weighted rotations of each
Lebedev grid.

Skala consumes `rho`, `grad(rho)`, and positive `tau`; it does not require a
density Hessian as an input tensor. Higher spatial derivatives nevertheless
enter the generalized Kohn-Sham operator through
`v_rho-div(dE/d grad(rho))` and `-0.5*div(v_tau*grad(psi))`. CP-PAW therefore
constructs these divergence operators with the same Fourier discretization
used for the corresponding primitive fields. This adjoint consistency is
required before the electronic operator can be used reliably in SCF or wave
function dynamics.

Atomic forces use the model's first reverse-mode derivatives together with
spatial derivatives of the PAW input fields. In particular, moving local grids
require `grad(rho)`, the density Hessian `grad(grad(rho))`, and `grad(tau)`.
These derivatives are generated by the native-grid interpolation and checked
against central finite differences. Derivatives of the model adjoints, i.e.
second model derivatives, are not needed for stationary-state forces; they
would be required for response properties such as phonons or force constants.
Unlike CP-PAW's nonspherical one-center XC Taylor expansion, which requests
second and third functional derivatives, the Skala path evaluates the complete
three-dimensional atom block directly on radial/Lebedev points.

Stationary-state stress follows the same distinction: it needs first model
derivatives and the spatial derivative structure above, including the kinetic
tensor associated with positive tau. Second derivatives of the Skala model are
not required unless a response derivative of stress or force is requested.

When `APPLY=T`, `SAFEORTHO` defaults to `F` because the robust
orthogonalization is required by the harder PAW one-center operator. An
explicit `SAFEORTHO` value is respected. Start electronic dynamics with a
conservative `DT**2/MPSI`; the conventional CP-PAW default is too aggressive
for the current experimental operator.

For PAW, the eventual caller must construct each atom block from primary fields
before inference:

```
smooth atom-partitioned field - pseudo one-center field + AE one-center field
```

"Atom-partitioned" means that each complete atom block receives its share of
the common quadrature rows and weights. The physical `rho`, `grad(rho)`, and
`tau` values on a retained row are not multiplied by the partition weight.
The local rows are assembled as `smooth + AE - pseudo` before the model call.

This differs from both separate conventional one-center XC corrections and a
literal copy of CP2K's GAPW implementation. Nonlinear Skala features are formed
only after the PAW fields have been combined.

Build and download the hash-pinned model with:

```sh
src/Buildtools/paw_skala_setup.sh --device auto --download-model
```

The setup script prints the installed bridge root and a smoke-test command.
The CP-PAW build can later consume that root through
`CPPAW_SKALA_FTORCH_ROOT` when the native PAW call path is enabled. The bridge
is composable with every existing library selection:

```sh
CPPAW_USE_SKALA_FTORCH=yes \
CPPAW_SKALA_FTORCH_ROOT="$PWD/bin/skala_ftorch_cpu" \
src/Buildtools/paw_build.sh -c fast -j16 -z
```

Use the CUDA bridge root with an NVIDIA target in the same way. Skala is never
enabled implicitly, so conventional CP-PAW binaries retain their dependency
set and behavior.

For CUDA, the setup script selects an `nvcc` whose major and minor toolkit
version matches the selected PyTorch package. Set `CUDACXX` to require a
specific compiler. CMake's CUDA architecture probe is initialized from
the selected GPU's compute capability; `CPPAW_SKALA_CUDA_ARCH` can override
that value on unusual or cross-compiled systems.

Binary PyTorch packages use the GNU OpenMP runtime, while `nvfortran` links the
NVIDIA OpenMP runtime even when CP-PAW does not enable OpenMP directives. The
NVHPC bridge therefore defaults Torch host work to one thread; override this
only with `CPPAW_SKALA_TORCH_THREADS`. Do not add `-mp` to a binary-PyTorch
build. A genuinely threaded NVHPC host configuration requires a LibTorch build
without GNU OpenMP rather than suppressing the mixed-runtime warning.

The optional periodic Si64 regression uses the conservative electronic
dynamics settings above and the standard Si64 structure:

```sh
cd tests/profile/si64
TEST=si64_skala SKALA_MODEL=/absolute/path/to/skala-1.1-rev1-cuda.fun \
  SKALA_DEVICE=AUTO SKALA_CHECK=F NSTEPS=1 CASES=gpu_all \
  ./run_benchmark.sh
```

Set `SKALA_CHECK=T` for variational diagnostics. The check mode is intended for
correctness runs; benchmark timings should use `SKALA_CHECK=F`.

`INTERPOLATEPARTITION=T` is an experimental performance mode. It interpolates
the local atom-grid weights from the cached exact smooth-grid partition. The
default remains `F` while energy and grid-convergence effects are evaluated.
The exact default path automatically reuses local partition weights for atoms
whose current periodic environments are related by a pure lattice translation
and whose augmentation cutoffs match.

`RADIALPOINTS`, `LEBEDEVEXACTNESS`, and `LEBEDEVORIENTATIONS` control the
moving local quadrature. Their defaults are 100, 17, and 1, respectively.
Additional orientations average deterministic rotations of the same Lebedev
rule and preserve the normalized quadrature weights. They expose and reduce
rotational integration error from nonlinear meta-GGA features. Energy,
particle number, forces, and stress should be converged with respect to all
three controls before production use; raising `RADIALPOINTS` is particularly
relevant for sharply peaked AE core fields.

The hybrid quadrature joins the moving radial/Lebedev PAW grid to the fixed
native cell grid with a quintic radial blend over the outer 20 percent of the
augmentation radius. Both the value and first derivative vanish at the
endpoints. Native-grid coordinates passed to Skala remain in the fixed cell
frame; minimum-image vectors are used only for the radial blend and PAW-local
quantities. This avoids a finite model-energy jump when a periodic atom crosses
an image-selection boundary.

Density, density-gradient, and kinetic-energy-density interpolation on each
local atom-grid point shares one native-grid stencil. The reverse mapping uses
the same weights as a combined discrete adjoint, which avoids recomputing the
stencil separately for the five primitive fields without changing the model
inputs or generalized-Kohn-Sham derivative. `CHECK=T` also reports energies,
row counts, and explicit force components per atom block. The complete force,
including the PAW projector response, remains experimental until end-to-end
finite differences agree at converged electronic states.
