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
`!CONTROL!DFT!SKALA` block activates primitive-field input preparation and the
positive-tau consistency check; it does not yet replace CP-PAW's conventional
XC energy. The protocol report labels this state explicitly.

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
