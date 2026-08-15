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
isolated with their component switches. Its safe default is `F`; forces and
stress are not yet implemented for this mode.

```text
!SKALA
 MODEL='path/to/skala-1.1-rev1-cuda.fun'
 DEVICE='AUTO'
 APPLY=F
 APPLYSMOOTH=T
 APPLYTAU=T
 APPLYONECENTER=T
 CHECK=F
!END
```

`CHECK=T` performs a one-time central finite-difference check of the model
adjoint and its PAW one-center density-matrix contraction. It also verifies
the discrete integration-by-parts identity for the density-gradient adjoint
and compares the occupied-state expectation of the positive-tau Hamiltonian
with the primitive `integral v_tau*tau`. `CHECK=F` omits the diagnostic field
copies and the additional wave-function overlap.

Skala consumes `rho`, `grad(rho)`, and positive `tau`; it does not require a
density Hessian as an input tensor. Higher spatial derivatives nevertheless
enter the generalized Kohn-Sham operator through
`v_rho-div(dE/d grad(rho))` and `-0.5*div(v_tau*grad(psi))`. CP-PAW therefore
constructs these divergence operators with the same Fourier discretization
used for the corresponding primitive fields. This adjoint consistency is
required before the electronic operator can be used reliably in SCF or wave
function dynamics.

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
PyTorch's supported architecture list; `CPPAW_SKALA_CUDA_ARCH` can override
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
