# Skala 1.1 FTorch bridge

This optional shared library isolates LibTorch and C++ from CP-PAW's legacy
Fortran build. It uses FTorch v1.0.0 for array-to-tensor transfers and a small
LibTorch protocol adapter for Skala's dictionary input and autograd output.
The setup path builds the pinned FTorch source even if another installation is
visible. CMake users can opt into an external package with
`-DCPPAW_SKALA_USE_SYSTEM_FTORCH=ON`.

### Linux BLAS compatibility

Skala-enabled Linux builds selected with `-lopenblas` link the static OpenBLAS
archive, including LAPACK, with `--exclude-libs,libopenblas.a`. The OpenBLAS
development package must supply this archive and its static dependencies via
`pkg-config`. This keeps CP-PAW's BLAS symbols local to its executable while
LibTorch uses its own BLAS. The ordinary non-Skala build is unchanged.

This is needed because LibTorch's exported BLAS symbols can collide with a
second implementation in the process; see
[PyTorch issue 182263](https://github.com/pytorch/pytorch/issues/182263).
Do not work around this by globally preloading `libtorch_cpu.so`: OpenBLAS
LAPACK and Torch's embedded complex-dot interface are not interchangeable.
The private archive avoids both directions of symbol preemption without
replacing CP-PAW's selected BLAS with Torch's embedded implementation.

### Functional contract

The scientific input contract is deliberately PAW-specific. Host arrays use
`density(npoint,2)`, `grad(npoint,3,2)`, and `kin(npoint,2)`. The bridge converts
them to Skala protocol-v2 tensors and returns derivatives in the host layout.
`kin` must be the positive kinetic-energy density
`0.5*sum_i |grad psi_i|^2`; CP-PAW's existing Laplacian-gauge `RHOKIN` is not a
drop-in replacement.

The native PAW call path is being integrated in stages. The
`!CONTROL!DFT!SKALA` block now assembles joint-source PAW atom blocks, evaluates
the model, and maps its density, density-gradient, and positive-tau adjoints
back to the one-center density matrices and smooth wave-function grid. The
smooth scalar operator includes the Fourier-space divergence of the gradient
adjoint. These operators are validated diagnostically; they do not replace
CP-PAW's conventional XC energy or Hamiltonian by default. `APPLY=T` enables
the experimental electronic operator, including the generalized Kohn-Sham
positive-tau term. The smooth scalar, positive-tau, and one-center parts can be
isolated with their component switches. Its safe default is `F`. Experimental
analytic atomic forces include the explicit model-coordinate, moving local
grid, interpolated-primitive, and atom-image partition terms in addition to
the PAW projector response. The smooth scalar adjoint is also contracted with
the translated pseudo-core densities after Skala inference; this contribution
cannot be taken from the earlier conventional-XC potential. Analytic stress is
available when CP-PAW requests cell stress. It includes the model-coordinate,
energy partition, independent descriptor-window, interpolated smooth-field,
all-source, pseudo-core, positive-tau, and existing PAW projector responses.
The new mapped radial grid has fixed Cartesian offsets and base weights, so
there is no affine local-quadrature volume factor or radial-blend derivative.

The functional is experimental. Kernel adjoint tests do not establish
electronic or quadrature convergence, nor end-to-end force/stress accuracy.

```text
!SKALA
 MODEL='path/to/skala-1.1-rev1.fun'
 DEVICE='CPU'
 RADIALPOINTS=200
 LEBEDEVEXACTNESS=53
 LEBEDEVORIENTATIONS=1
 IMAGESHELLS=1
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

`CHECK=T` reports the one-center density-matrix/primitive contraction identity.
The source and image-kernel unit tests separately check their finite-difference
derivatives. The applied-functional check also verifies
the discrete integration-by-parts identity for the density-gradient adjoint
and compares the occupied-state expectation of the positive-tau Hamiltonian
with the primitive `integral v_tau*tau`. `CHECK=F` omits the diagnostic field
copies and the additional wave-function overlap. Force diagnostics report the
total energy at full double precision and list the pseudo-core contribution
separately, so central differences can be evaluated without the rounded energy
summary. Stress diagnostics report each Skala contribution and the complete
`D E / D STRAIN` tensor before CP-PAW converts it to its internal stress sign
convention.

Check mode also reports the PAW-metric electronic residual
`R_j = H psi_j - S sum_i psi_i <psi_i|H|psi_j>`, its occupied RMS and maximum,
the occupation commutator `|H_ij| |f_i-f_j| / max(f)`, `max|Psi^dagger S Psi-I|`,
and Hamiltonian Hermiticity. These are diagnostics of the current Hamiltonian:
with `APPLY=F` they refer to conventional XC, not Skala. The residual assumes
an orthonormal PAW basis; its separate overlap check must pass. Small energy
changes or a stationary subspace alone are not a sufficient convergence test.
Even passing all these checks does not prove the global ground-state minimum.
The [stationarity checker](../../tests/fulltests/skala_si2/stationarity.py)
checks every step for finite values and metric/operator consistency; optional
explicit residual and occupation-commutator tolerances test the final window.

The one-center tau audit separately compares the angular positive-tau integral
with an angularly exact radial contraction of the same partial waves, density
matrix and radial domain. This is valence tau; production reconstruction adds
the frozen-core contribution separately. Further diagnostics expose the outer
radial domain, surface term, nominal ZORA contribution and a numerical radial
kinetic operator. The legacy HBS setup matrix is constructed using its stored
equation-based kinetic action and replaced pseudo-wave tails, so the raw
`TAU-KINETIC DIFFERENCE` is not an identity test. No correction or rescaling of
the Skala tau input is made to force it to match that matrix.

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
than an affine deformation of the radial grid. The Si2 test drivers cover
isotropic, uniaxial, and symmetric-shear strains and MPI parity. Use a state
converged with the same functional and reconstruction as the test executable.

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
For CPU execution of a CUDA-capable binary, set both `CPPAW_GPU_MODE=off`
for the PAW acceleration libraries and `DEVICE='CPU'` in `!SKALA` (or
`SKALA_DEVICE=CPU` in these test drivers), using the CPU model export.
The two controls are independent: disabling PAW GPU acceleration alone does
not prevent an `AUTO` Skala model from selecting CUDA. This runtime CPU mode
is distinct from the GPU-library-free GNU installation described below.
`SKALA_MPI_NORM_TOLERANCE` controls the separate absolute tolerance for the
large smooth-operator L2 diagnostics and defaults to `2e-5`.
`quadrature_convergence.sh` defaults to radial counts `100 200 400` at Lebedev
exactness 17 to expose coarse-grid errors and reports every result relative to
the final, finest case. The
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

For PAW, the caller constructs each target atom block from primary fields
before inference:

```
smooth field + sum_(source atoms and images) (AE source - pseudo source)
```

The atom label identifies the target descriptor block, not the only source of
its reconstructed fields. Every source within its radial support contributes
to every target row, including periodic images. The smooth density contains
pseudo core; the source density adds AE minus pseudo core. Smooth positive tau
contains valence only, so the complete frozen-core positive tau is added once.
It is derived from the occupied setup core orbitals, not the Laplacian-gauge
kinetic-energy density. Collinear total/spin fields are converted to up/down
channels after reconstruction; the core is unpolarized.

Energy weights are `w*p_A0`, where `p_A0` uses the target-specific fixed image
layout of all atoms. Descriptor weights are `w*T(d_A)`, with `d_A` recomputed
from the target's self images only and a quintic taper from `1e-12` to `1e-11`.
The descriptor window must not be obtained by tapering the all-atom energy
partition. Rows with zero energy weight can still contribute to descriptors.
Neither weight multiplies the physical rho, gradient or tau. The partition
uses shifts +/-1 around each source image nearest the target; field images
are enumerated separately from the actual source support, not that fixed shell.

All source density-matrix adjoints are accumulated across target blocks and
MPI ranks before building any one-center Hamiltonian. The reconstruction,
matrix adjoint and spatial derivatives use the same radial interpolation
polynomial. Forces include target and source motion; strain also includes
source-image translations and the independent self-image window derivative.

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
The dedicated targets retain CP-PAW's host FFT backend: bundled NVPL FFTW is
used where present, while x86 NVHPC installations without NVPL require a
normal FFTW3 development installation discoverable through `pkg-config`.
NVHPC otherwise uses the CUDA toolkit bundled with the compiler. On a host
whose driver supports only an older co-installed toolkit, set
`NVHPC_CUDA_HOME` before a clean CP-PAW build and build the CUDA FTorch bridge
for the same CUDA major/minor version. This keeps both OpenACC device images
and LibTorch below the driver's supported CUDA level.

Set `CPPAW_SKALA_GRID_BACK_ACC=1` to use it; the default is off. The output
density, gradient, and kinetic-energy-density adjoint grids remain resident
from the forward-begin boundary through all atom blocks and return to the host
at forward-end. Smooth points use a conflict-free device loop, while local
eighth-order interpolation stencils use atomic accumulation. Set
`CPPAW_SKALA_GRID_BACK_ACC_MIN_POINTS` to override the default threshold of
32768 smooth-grid cells. The CPU batch path remains available in every build,
including builds without OpenACC or CUDA. Profile builds report kernel time as
`SKALA_GRID_BACK_ACC_SCATTER` and transfers as
`ACC_COPY_SKALA_GRID_BACK_{IN,INPUT,OUT}`.
In the default MPI-root CUDA mode, only the root rank creates this residency;
distributed CPU or CUDA inference enables it only on ranks that own atom
blocks.

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
One host thread does not itself resolve the runtime conflict or certify its
safety; see [NVIDIA's OpenMP compatibility warning](https://docs.nvidia.com/nvpl/latest/).

Set `CPPAW_SKALA_DETERMINISTIC=1` to request PyTorch's deterministic algorithm
mode, disable TF32, and install the deterministic cuBLAS workspace setting.
This is opt-in and does not guarantee identical model adjoints across processes
or hardware. MPI-root inference remains the default.

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

`INTERPOLATEPARTITION=T` is rejected by the joint-source implementation. The
old smooth-grid partition cache cannot represent its two independent weight
families. `F` is accepted for input compatibility. Exact weights are cached for
an unchanged geometry and cell, with no source-field caching.

`RADIALPOINTS`, `LEBEDEVEXACTNESS`, and `LEBEDEVORIENTATIONS` control the
moving local quadrature. Their starting defaults are 200, 53, and
1, respectively, not a convergence guarantee. Additional orientations average
deterministic rotations of the same Lebedev rule and preserve the normalized
quadrature weights. They expose and reduce rotational integration error from
nonlinear meta-GGA features.

Angular exactness is a minimum: the library selects the next supported rule,
up to 65. The largest rules are 53 (974 points), 59 (1202 points), and 65
(1454 points). In particular, `LEBEDEVEXACTNESS=64` selects order 65; the
protocol reports both the requested and actual exactness and the angular
point count. At 200 radial points and one orientation this creates 290800
candidate rows per atom, versus 194800 for order 53, before zero-weight
filtering. The roughly 49% increase in quadrature rows also increases model
work and memory demand. It does not replace native-grid cutoff convergence.
See [`LEBEDEV.md`](../../tests/unittests/skala_reconstruction/LEBEDEV.md) for
the coefficient provenance and polynomial tests.

Skala's nonlocal atom layers require complete atom blocks, so peak accelerator
memory grows with the number of orientations. Use 200/53/1 as a starting point
and converge energy, particle number, forces, and stress explicitly for
demanding calculations.

The corrected quadrature uses Gauss-Legendre nodes mapped to `[0,infinity)`
with `r=x/(1-x)` in bohr and Lebedev angular rules. Native smooth fields are
interpolated onto these rows. There is no owner-sphere cutoff on foreign source
fields, radial blend or direct native-cell energy quadrature. Both quadrature
and native-grid interpolation must be converged, especially around foreign
nuclei. A finite electron-count error is reported rather than normalized away.

`IMAGESHELLS` independently controls the compact periodic support for the
energy partition and self-image descriptor window. Its default of one shell
is a starting point, not a converged periodic quadrature. The reconstructed
source fields instead include every image within their radial support.
Converge the image support independently of radial/angular resolution.
The protocol reports the constant-field integral (`PARTITION VOLUME`), the
exact cell volume, and their relative difference, with a warning above `1e-4`.
A smaller volume error alone does not certify energy/force convergence. No
density or weight normalization is applied to hide an integration error.

Use the model-free periodic probe to check the quadrature independently of
Skala inference and PAW density reconstruction:

```sh
make -C bin/Build_fast -f Makefile \
  -f "$PWD/tests/unittests/skala_reconstruction/driver.mk" skala-partition-probe
# Arguments: image shells, radial points, angular exactness, orientations,
# optional relative tolerance (nonzero exit if volume or first mode fails).
bin/Build_fast/unit-tests/partition_measure.x 1 200 53 1 1e-4
```

Run the reconstruction kernels after a normal build (no model is required):

```sh
bash tests/unittests/skala_reconstruction/run.sh bin/Build_fast
```

The [crystal probes](../../tests/fulltests/skala_crystals/README.md) use the
CO2, NH3 and urea structures from the CP2K manuscript benchmark repository.
They test CP-PAW execution and adjoints, not agreement between PAW and GAPW
total-energy zeros.

Density, density-gradient, and kinetic-energy-density interpolation on each
local atom-grid point shares one native-grid stencil. The reverse mapping uses
the same weights as a combined discrete adjoint, which avoids recomputing the
stencil separately for the five primitive fields without changing the model
inputs or generalized-Kohn-Sham derivative. `CHECK=T` also reports energies,
row counts, and explicit force components per atom block. The complete force,
including the PAW projector response, remains experimental until end-to-end
finite differences agree at converged electronic states.
