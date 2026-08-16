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
 MODEL='path/to/skala-1.1-rev1.fun'
 DEVICE='CPU'
 RADIALPOINTS=200
 LEBEDEVEXACTNESS=53
 LEBEDEVORIENTATIONS=1
 APPLY=F
 APPLYSMOOTH=T
 APPLYTAU=T
 APPLYONECENTER=T
 INTERPOLATEPARTITION=F
 DISTRIBUTEGPU=F
 CHECK=F
!END
```

CUDA model inference defaults to MPI-root execution. This keeps complete PAW
atom blocks in one LibTorch instance and makes energy, force, and stress
independent of how MPI distributes atoms. `DISTRIBUTEGPU=T` restores the
experimental per-rank CUDA path. It can expose more GPU parallelism, but the
Skala 1.1 float32 network produced rank-dependent adjoints when several MPI
processes evaluated different atoms, so it is not the correctness default.
CPU model inference remains distributed across MPI ranks.

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
`SKALA_DEVICE=CPU|CUDA|AUTO` selects the model device without hiding CUDA from
an otherwise GPU-enabled CP-PAW executable.
`SKALA_MPI_NORM_TOLERANCE` controls the separate absolute tolerance for the
large smooth-operator L2 diagnostics and defaults to `2e-5`.
`quadrature_convergence.sh` defaults to radial counts `100 200 400` at Lebedev
exactness 17 to expose coarse-grid errors and reports every result relative to
the final, finest case. The
lists can be changed with `SKALA_RADIAL_POINT_LIST` and
`SKALA_LEBEDEV_EXACTNESS_LIST`. `SKALA_LEBEDEV_ORIENTATION_LIST` adds a
convergence sweep over deterministic, equally weighted rotations of each
Lebedev grid.

The 2026-08-15 end-to-end validation used a stationary Si2 restart. At
`200/53/1`, the selected force component differed from its central finite
difference by `8.55e-4 H/bohr`. Isotropic, uniaxial, and symmetric-shear
stress checks differed by `7.39e-4`, `8.78e-5`, and `4.95e-4 H`, respectively;
all were below the `2e-3` absolute tolerance. One- and four-rank Gamma runs
agreed within `4.1e-13 H` in energy, `3.8e-8 H/bohr` in forces, and
`1.7e-7 H` in stress on the tested CPU/GPU paths. The periodic eight-k-point
integration test also passed. These checks exercise total energies and
stationary-state derivatives, not just the isolated Torch protocol.

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

For a CPU-only installation, first provide a CPU PyTorch or LibTorch package,
then build the bridge and download the hash-pinned CPU model with GNU Fortran:

```sh
FC=gfortran src/Buildtools/paw_skala_setup.sh --device cpu --download-model
```

The setup script prints the installed bridge root and a smoke-test command.
It fetches the pinned FTorch source, but never installs PyTorch or LibTorch.
The dedicated CPU targets select the GNU toolchain and the default
`bin/skala_ftorch_cpu` root automatically:

```sh
src/Buildtools/paw_build.sh -c skala_cpu_fast -j16 -z
src/Buildtools/paw_build.sh -c skala_cpu_fast_parallel -j16 -z

# Or build both through the installer without probing NVIDIA targets.
CPPAW_INSTALL_NVHPC=no CPPAW_INSTALL_SKALA_CPU=require ./paw_install
```

Set `CPPAW_SKALA_FTORCH_ROOT` for a custom CPU bridge prefix. Use a CUDA bridge
root with an NVIDIA target in the same way. `CPPAW_INSTALL_SKALA_CPU` defaults
to `no`, and the conventional `dbg`, `fast`, and `fast_parallel` binaries never
enable Skala implicitly, so they retain their dependency set and behavior.

The dedicated `nvhpc_skala_grid_acc_{fast,profile}` targets add the optional
OpenACC native-grid adjoint kernel. They still require
`CPPAW_USE_SKALA_FTORCH=yes` and a CUDA FTorch root for model inference. The
combined `nvhpc_gpu_acc_*` and `nvhpc_gpu_all_*` targets already compile this
kernel because they enable OpenACC for the other accelerator paths.

Set `CPPAW_SKALA_GRID_BACK_ACC=1` to use it; the default is off. The output
density, gradient, and kinetic-energy-density adjoint grids remain resident
from the forward-begin boundary through all atom blocks and return to the host
at forward-end. Smooth points use a conflict-free device loop, while local
eighth-order interpolation stencils use atomic accumulation. Set
`CPPAW_SKALA_GRID_BACK_ACC_MIN_POINTS` to override the default threshold of
32768 smooth-grid cells. The CPU batch path remains available in every build,
including builds without OpenACC or CUDA. Profile builds report kernel time as
`SKALA_GRID_BACK_ACC_SCATTER` and transfers as
`ACC_COPY_SKALA_GRID_BACK_{RESIDENCY_IN,INPUT,RESIDENCY_OUT}`.

For CUDA, the setup script selects an `nvcc` whose major and minor toolkit
version matches the selected PyTorch package. Set `CUDACXX` to require a
specific compiler. It probes the CUDA host C++ compiler and, when necessary,
selects an installed compatible GNU version. `CUDAHOSTCXX` makes that choice
explicit. A narrowly scoped `-U_GNU_SOURCE` fallback handles CUDA 12.4 with
glibc 2.41 headers after a compile probe confirms that combination. CMake's
CUDA architecture probe is initialized from the selected GPU's compute
capability; `CPPAW_SKALA_CUDA_ARCH` can override that value on unusual or
cross-compiled systems. Custom installation prefixes receive independent
build trees; `CPPAW_SKALA_BUILD_DIR` can select one explicitly.

Binary PyTorch packages use the GNU OpenMP runtime, while `nvfortran` links the
NVIDIA OpenMP runtime even when CP-PAW does not enable OpenMP directives. The
NVHPC bridge therefore defaults Torch host work to one thread; override this
only with `CPPAW_SKALA_TORCH_THREADS`. Do not add `-mp` to a binary-PyTorch
build. A genuinely threaded NVHPC host configuration requires a LibTorch build
without GNU OpenMP rather than suppressing the mixed-runtime warning.

Set `CPPAW_SKALA_DETERMINISTIC=1` to request PyTorch's deterministic algorithm
mode, disable TF32, and install the deterministic cuBLAS workspace setting.
This is opt-in because it did not remove cross-process model-adjoint variation
for Skala 1.1 and slowed the tested CUDA path. MPI-root inference is the
reproducible default instead.

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
moving local quadrature. Their portable production defaults are 200, 53, and
1, respectively. A Si2 sweep found the former 100/17/1 grid too coarse. At
eight orientations, 300, 360, and 400 radial points agreed within about
`5e-6` Hartree in model XC energy. Additional orientations average
deterministic rotations of the same Lebedev rule and preserve the normalized
quadrature weights. They expose and reduce rotational integration error from
nonlinear meta-GGA features.

The 300/53/8 setting is a memory-intensive reference grid: Skala's nonlocal
atom layers require complete atom blocks, so peak accelerator memory grows
with the number of orientations. It required about 74 GB on the tested CUDA
build. Use 200/53/1 as the portable starting point and converge energy,
particle number, forces, and stress explicitly for demanding calculations.

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
