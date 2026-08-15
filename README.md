<p align="center">
<a href="https://cppaw.org">
<img src="src/Docs/Figs/PAWlogo/paw_github.svg" width="300" title="cppaw.org">
</a>
</p>
<p align="right"> 
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="License: GPLv3"></a>
</p>

# CP-PAW code package

See https://cppaw.org for further information (Currently, the
description on https://cppaw.org refers to an older release and does
not apply to the present implementation.)

## `cp-paw-nvhpc` development features

> [!WARNING]
> The NVIDIA accelerator paths and the Skala functional on this branch are
> experimental. CPU fallbacks remain available, and optional accelerator
> features are either selected by a dedicated build target or guarded by a
> run-time switch. Validate energies, forces, stress, and parallel parity for
> the intended system before using these paths for production calculations.

The build system detects the NVIDIA HPC SDK, CUDA, and individual libraries
instead of assuming that every NVIDIA installation provides the same stack.
With the installer's default `auto` policy, unavailable optional targets are
skipped. Set an installer option to `require` to turn a missing dependency or
failed build into an installation error.

### NVIDIA library coverage

| Component | CP-PAW integration | Build targets |
| --- | --- | --- |
| NVPL | CPU BLAS, LAPACK, and FFTW backend when found; the NVHPC compiler BLAS/LAPACK libraries remain a fallback | `nvhpc_fast*`, `nvhpc_profile*` |
| NVBLAS | BLAS-3 interposition experiment for existing DGEMM, ZGEMM, DSYRK, and ZHERK calls | `nvhpc_nvblas_*` |
| NVLAMATH | NVIDIA LAPACK/cuSOLVER wrapper path | `nvhpc_nvlamath_*`, `nvhpc_gpu_all_*` |
| cuFFTW | FFTW3-compatible cuFFT wrapper that preserves the existing CP-PAW FFT call structure | `nvhpc_cufftw_*`, `nvhpc_gpu_all_*` |
| cuFFT | Native cuFFT/OpenACC path for selected batched complex FFTs; the 1-D and diagnostic 3-D paths are opt-in at run time | `nvhpc_cufft_*`, `nvhpc_gpu_acc_*`, `nvhpc_gpu_all_*` |
| cuBLAS | Explicit OpenACC/cuBLAS path for selected dense matrix, overlap, orthogonalization, projection, and one-center operations | `nvhpc_cublas_acc_*`, `nvhpc_gpu_acc_*`, `nvhpc_gpu_all_*` |
| cuSOLVER | Standard and generalized real/complex eigensolvers plus an opt-in Gram-Cholesky path | `nvhpc_cusolver_acc_*`, `nvhpc_gpu_acc_*`, `nvhpc_gpu_all_*` |
| FTorch/LibTorch | Optional bridge to the experimental Skala 1.1 PAW functional | composable with every build target |

The combined `nvhpc_gpu_acc_*` targets enable the native cuFFT, explicit
cuBLAS, and explicit cuSOLVER integrations in one executable. The optional
`nvhpc_gpu_all_*` targets additionally link cuFFTW and NVLAMATH so that the
complete implemented library stack can be compiled and tested together.
NVBLAS remains separate because it interposes the host BLAS interface.

```sh
# Combined release builds
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_fast -j16 -z
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_fast_parallel -j16 -z

# Combined profiling build with the device-residency defaults
CPPAW_TOOLCHAIN=nvhpc \
  src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_residency_profile -j16 -z

# Integration build containing all composable NVIDIA libraries
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_fast -j16 -z
```

Important run-time defaults are deliberately conservative:

- explicit cuBLAS is enabled in builds that contain it, with an offload
  threshold of `CPPAW_CUBLAS_ACC_MINFLOP=1e7`;
- cuSOLVER is enabled for supported problems of size 256 or larger and can be
  disabled with `CPPAW_CUSOLVER_ACC=0`;
- native cuFFT is disabled until `CPPAW_CUFFT_ACC=1` is set; its default 1-D
  threshold is 512 batched elements, while the native 3-D path remains a
  separate diagnostic option;
- focused residency defaults are enabled by the
  `nvhpc_gpu_acc_residency_profile*` targets; the broader all-in-one stack is
  opt-in via `CPPAW_GPU_RESIDENCY_STACK=1`;
- Skala application and cuBLAS FP64 emulation are disabled by default.

The all-library builds do not add cuFFTMp, cuBLASMp, cuSOLVERMp, cuEST, or
ALCHEMI. In particular, the Ozaki implementation described below uses the
CUDA 13 cuBLAS API and has no cuEST dependency.

### CUDA 13 FP64 Ozaki emulation

CUDA 13.0 Update 2 introduced a cuBLAS fixed-point emulation API for FP64
matrix operations. CP-PAW detects this API with a compile probe and, when it is
available, compiles an opt-in dynamic-mantissa Ozaki path into every explicit
cuBLAS target. Older CUDA toolkits and unsupported operations continue to use
native FP64 without another source or link dependency.

```sh
export CPPAW_CUBLAS_FP64_EMULATION=1
export CPPAW_CUBLAS_FP64_STRATEGY=performant
export CPPAW_CUBLAS_FP64_WORKSPACE_MB=2048
```

Safe per-kernel defaults after the main switch is enabled are
`CPPAW_CUBLAS_FP64_DGEMM=1`, `CPPAW_CUBLAS_FP64_ZGEMM=1`, and
`CPPAW_CUBLAS_FP64_ZHERK=0`. Dynamic mantissa control, the `performant`
dispatch strategy, and one persistent 2048 MiB cuBLAS workspace are used. The
workspace is rebound after stream changes. `CPPAW_CUBLAS_FP64_STRATEGY=eager`
is intended for capability and correctness checks rather than routine use.

Set `CPPAW_CUBLAS_FP64_TELEMETRY=1` to record whether each selected operation
used emulation, fell back to native FP64, or could not report its mode. Set
`CPPAW_CUBLAS_FP64_CHECK=1` for machine-readable energy, orthonormality, force,
stress, and k-point validation data. These diagnostics synchronize additional
GPU work and should be disabled for timing runs. The complete switch reference
and the 1024/2048/4096-band validation driver are documented in
[`tests/profile/README.md`](tests/profile/README.md) and implemented by
[`tests/profile/si64/run_ozaki_bands.sh`](tests/profile/si64/run_ozaki_bands.sh).

### Experimental Skala 1.1 PAW functional

The optional Skala path uses FTorch and LibTorch to evaluate a TorchScript
Skala 1.1 model on CP-PAW's native smooth grid and atom-centered PAW grids. It
assembles each atom block as smooth plus all-electron one-center minus pseudo
one-center fields before model inference. This is a PAW-specific integration;
it is not a literal copy of CP2K's GAPW implementation.

Skala consumes the density, its gradient, and the positive kinetic-energy
density. CP-PAW constructs the corresponding generalized Kohn-Sham scalar and
positive-tau operators, one-center contributions, and the spatial derivatives
needed for analytic forces and stress. Higher spatial derivatives enter those
operator and moving-grid contractions even though the model input itself does
not contain a density Hessian.

Build a pinned FTorch bridge and download the model, then compose the bridge
with the desired CP-PAW target:

```sh
src/Buildtools/paw_skala_setup.sh --device auto --download-model

CPPAW_USE_SKALA_FTORCH=yes \
CPPAW_SKALA_FTORCH_ROOT="$PWD/bin/skala_ftorch_cuda" \
  src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_fast -j16 -z
```

The functional is selected in the input with a `!SKALA` block. The conservative
starting point evaluates the Skala path while keeping CP-PAW's conventional XC
operator active:

```text
!SKALA
 MODEL='path/to/skala-1.1-rev1-cuda.fun'
 DEVICE='AUTO'
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

`APPLY=T` enables the experimental electronic operator. CUDA inference runs on
the MPI root by default because that is the reproducible multi-rank mode for
the current float32 model; `DISTRIBUTEGPU=T` is experimental. `CHECK=T` enables
variational diagnostics and should be used for correctness runs. Complete
setup notes, SCF cautions, quadrature guidance, force/stress finite-difference
drivers, MPI and k-point validation, and the exact field contract are in
[`src/SkalaBridge/README.md`](src/SkalaBridge/README.md).

### Profiling and validation

Profiling builds record FFT, BLAS/LAPACK, accelerator, MPI transpose, transfer,
residency, and Skala phase data in `cppaw_accel_profile*.csv`. The Si64 harness
can compare CPU and GPU configurations and emit machine-readable TSV plus
Markdown summaries. These cases are intentionally outside the default test
suite and make no general performance guarantee.

```sh
src/Tools/Scripts/paw_gpu_capabilities.sh
cd tests/profile/si64
NSTEPS=1 ./run_gpu_exploration.sh
```

See [`tests/profile/README.md`](tests/profile/README.md) for the build/case
matrix, run-time switches, correctness checks, and interpretation of profiler
columns.


# Configuration and Installation instructions


> [!IMPORTANT]
> The installation instructions are for the current version of the CP-PAW code. The instructions may not be applicable to older versions of the code. Installation process is under current development. Please report any issues to the developers.

## Requirements

- fortran compiler (e.g. gfortran, ifort; NVIDIA HPC SDK/nvfortran is supported)
- pkg-config, latexmk, GNU Make (version 4.3 or later)
- bash, cpp, ar
- tex (latex) distribution (e.g. TeX Live) 
- LAPACK, BLAS, FFTW3, MPI (optional), LIBXC (optional)
- optional NVIDIA HPC SDK and CUDA libraries, as summarized in the
  [`cp-paw-nvhpc` development features](#cp-paw-nvhpc-development-features)
  section above
- optional CMake, Python, FTorch, and LibTorch for the Skala bridge; the setup
  helper installs the pinned bridge dependencies
- tools: xmgrace, gnuplot, avogadro1

## Installation

1. Download the cppaw distribution from [Github](https://github.com/cp-paw/cp-paw).
2. Unpack the distribution in a directory. I will refer to this directory as the base directory.
3. Add the following lines to your profile (e.g. ~/.zshrc, ~/.bashrc, ~/.profile). Replace the definition of `PAWDIR` by the name of the base directory.
   ```
   export PAWDIR="path to your cppaw distribution"
   export PATH=${PAWDIR}/bin/fast:${PAWDIR}/bin/fast_parallel:${PAWDIR}/bin/dbg:${PATH}
   ```
4. Execute the following commands in the base directory:
   ```
   ./paw_install
   ```
   On systems where the NVIDIA HPC SDK is available, the installer will additionally try CPU/NVPL builds. Use:
   ```
   CPPAW_INSTALL_NVHPC=require ./paw_install
   ```
   to make those builds mandatory. CUDA-dependent NVBLAS, NVLAMATH, cuFFTW, native cuFFT, combined GPU, all-library GPU, cuBLAS/OpenACC and cuSOLVER/OpenACC variants are only attempted when CUDA is detected, unless requested explicitly with `CPPAW_INSTALL_NVBLAS=require`, `CPPAW_INSTALL_NVLAMATH=require`, `CPPAW_INSTALL_CUFFTW=require`, `CPPAW_INSTALL_CUFFT=require`, `CPPAW_INSTALL_GPU_ACC=require`, `CPPAW_INSTALL_GPU_ALL=require`, `CPPAW_INSTALL_CUBLAS_ACC=require` or `CPPAW_INSTALL_CUSOLVER_ACC=require`.
   With `CPPAW_INSTALL_PROFILE=yes`, the installer now also attempts the
   recommended GPU residency profile when CUDA/cuFFT/cuBLAS/cuSOLVER are
   available. Set `CPPAW_INSTALL_GPU_RESIDENCY_PROFILE=no` to skip it, or
   `CPPAW_INSTALL_GPU_RESIDENCY_PROFILE=require` to make it mandatory.
   The NVIDIA builds can also be selected directly:
   ```
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvlamath_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvlamath_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_cublas_acc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_cublas_acc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_fast_parallel
   ```
   Profiling builds add CP-PAW hotspot instrumentation for FFT, BLAS/LAPACK-style kernels and MPI transposes:
   ```
   src/Buildtools/paw_build.sh -c profile
   src/Buildtools/paw_build.sh -c profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvlamath_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvlamath_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_cublas_acc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_cublas_acc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_residency_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_residency_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_managed_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_unified_profile
   ```
   Profiled runs write `cppaw_accel_profile.csv` in serial mode and `cppaw_accel_profile.rankNNNNN.csv` in parallel mode. Set `CPPAW_ACCEL_PROFILE_FILE=<prefix>` to choose another file prefix, or `CPPAW_ACCEL_PROFILE=0` to disable collection at run time. In `nvhpc_cufft_*` and `nvhpc_gpu_acc_*` builds, the native cuFFT path is opt-in: set `CPPAW_CUFFT_ACC=1` to enable it; the profiled default offload threshold is 512 batched elements and `CPPAW_CUFFT_ACC_MIN_ELEMENTS=<elements>` overrides it. `CPPAW_CUFFT_ACC_3D=1` additionally enables the experimental 3-D cuFFT wrapper for diagnostics. Set `CPPAW_FFT_SERIAL_3D=1` only for single-rank diagnostics that route `PLANEWAVE$FFT` through a full-grid 3-D transform; the benchmark case is `gpu_resident_stack_serial3dfft`. `CPPAW_FFT_SERIAL_3D_ACC_MAP=1` additionally tests device-side sparse/full-grid mapping and is intentionally a separate diagnostic case (`gpu_resident_stack_serial3dfft_accmap`) because it is architecture-sensitive. Combined GPU targets use explicit NVHPC GPU memory mode by default (`CPPAW_NVHPC_GPU_MEMORY_MODE=separate`); the `nvhpc_gpu_acc_managed_profile` and `nvhpc_gpu_acc_unified_profile` targets build `-gpu=mem:managed` and `-gpu=mem:unified` variants for device-residency experiments. The `nvhpc_gpu_acc_residency_profile*` targets default to `CPPAW_GPU_RESIDENCY=1`; set `CPPAW_GPU_RESIDENCY=0` to run the same binary with residency disabled. Residency lets cuBLAS scalarproducts and addproducts reuse OpenACC-present operands where wavefunction loops already keep data on the device; eligible projection loops, including superwave cases with the gamma correction applied on device, also scatter `PROPSI` on the device, expand per-atom projector blocks on the GPU, and copy the final projection array back once, the setup/orthogonalization path can keep `PSI0` and phase-eligible `PSIM` resident across Gram-Schmidt, projection, propagation, and pseudo-overlap boundaries (`CPPAW_GPU_SETUP_PSIM_RESIDENCY=0` disables only the setup `PSIM` extension), and the opt-in HPSI diagnostic keeps ETOT-local `HPSI`/`PSI0` data present for the immediate expectation/Hamiltonian overlaps. Set `CPPAW_GPU_PRO_EXPANSION=0` to keep projector expansion on the host for diagnostics. In `nvhpc_cublas_acc_*` and `nvhpc_gpu_acc_*` builds, `CPPAW_CUBLAS_ACC=0` disables the cuBLAS path at run time and `CPPAW_CUBLAS_ACC_MINFLOP=<flops>` adjusts the offload threshold. The opt-in `CPPAW_GPU_DENMAT_ENERGY=1` diagnostic offloads the time-inversion one-center DENMAT energy/Lambda contraction; `CPPAW_GPU_DENMAT_MINFLOP=<flops>` adjusts that threshold. In `nvhpc_cusolver_acc_*` and `nvhpc_gpu_acc_*` builds, `CPPAW_CUSOLVER_ACC=0` disables the cuSOLVER path at run time and `CPPAW_CUSOLVER_ACC_MIN_N=<n>` adjusts the dense eigensolver offload threshold; `CPPAW_CUSOLVER_MIN_N=<n>` is accepted as a shorter alias. The default cuBLAS threshold is `1e7`; the Si64 profile favors this over more conservative thresholds. The default cuSOLVER threshold is `256`; smaller Si64 profiling runs can set `CPPAW_CUSOLVER_ACC_MIN_N=1` to force offload. Set `CPPAW_CUSOLVER_ACC_CHECK=1` to validate cuSOLVER eigensolver results by residual and orthonormality before accepting them; failures fall back to the CPU LAPACK path.
   To include profile targets in the installer run:
   ```
   CPPAW_INSTALL_PROFILE=yes ./paw_install
   ```
   The installer compiles with `CPPAW_INSTALL_JOBS=16` by default; set another
   value if your build host needs a smaller or larger parallel make. Set
   `CPPAW_INSTALL_GPU_ALL=yes` to add the optional all-library GPU binaries, or
   `CPPAW_INSTALL_GPU_MEMORY_PROFILES=yes` to add the optional managed/unified
   GPU profile binaries.
   A reusable 64-atom periodic silicon profiling case is available under
   `tests/profile/si64`. It is intentionally not part of the default test
   suite:
   ```
   cd tests/profile/si64
   PAWX="mpirun -np 4 ../../../bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" make all
   ```
   To try NVIDIA GPU BLAS interposition for the existing BLAS-3 calls, build an
   `nvhpc_nvblas_*` target and run the same case with `NVBLAS=yes`. The helper
   script `paw_nvblas.sh` creates an NVBLAS configuration file when
   `NVBLAS_CONFIG_FILE` is not already set:
   ```
   cd tests/profile/si64
   PATH="../../../bin/nvhpc_nvblas_profile_parallel:${PATH}" \
   PAWX="mpirun -np 4 ../../../bin/nvhpc_nvblas_profile_parallel/ppaw_nvhpc_nvblas_profile.x" \
   NVBLAS=yes make all
   ```
   To try NVLAMATH's LAPACK/cuSOLVER wrapper path for the existing LAPACK calls:
   ```
   cd tests/profile/si64
   PAWX="../../../bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" \
   make all
   ```
   To try cuFFTW/cuFFT for CP-PAW's existing FFTW3 calls:
   ```
   cd tests/profile/si64
   PAWX="mpirun -np 4 ../../../bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" \
   make all
   ```
   To try the native cuFFT/OpenACC 1-D FFT path:
   ```
   cd tests/profile/si64
   CPPAW_CUFFT_ACC=1 \
   PAWX="mpirun -np 4 ../../../bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" \
   make all
   ```
   For explicit cuBLAS/OpenACC profiling without NVBLAS interposition:
   ```
   cd tests/profile/si64
   PAWX="mpirun -np 4 ../../../bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" \
   make all
   ```
   For the recommended combined GPU residency profiling path:
   ```
   cd tests/profile/si64
   PAWX="../../../bin/nvhpc_gpu_acc_residency_profile/paw_nvhpc_gpu_acc_residency_profile.x" \
   make all
   ```
   To scan GPU/NVIDIA-library capabilities and run a short diagnostic matrix:
   ```
   src/Tools/Scripts/paw_gpu_capabilities.sh
   cd tests/profile/si64
   NSTEPS=1 ./run_gpu_exploration.sh
   ```
   The capability helper reports host FFTW library/include and BLAS/LAPACK
   availability plus current `recommended_*` case lists. On CUDA systems with
   cuBLAS and a usable host numerical stack it points routine GPU checks at the
   residency stack rather than the older force-all diagnostic cases. The
   standard, exploration, and overnight benchmark wrappers consume those
   recommendations through their `auto` case-list defaults.
   To force all native paths for diagnostics, add `CPPAW_CUFFT_ACC=1` and
   `CPPAW_CUSOLVER_ACC_MIN_N=1`.
   To inspect a benchmark matrix without starting CP-PAW runs, use the dry-run
   mode. It writes `metadata.txt` with the executable, MPI launcher, runtime
   environment, accelerator case variables, and planned command line for each
   selected case:
   ```
   cd tests/profile/si64
   DRY_RUN=yes CASES="gpu_resident_stack gpu_no_cufft" ./run_benchmark.sh
   ```
   For explicit cuSOLVER/OpenACC dense eigensolver profiling:
   ```
   cd tests/profile/si64
   CPPAW_CUSOLVER_ACC_MIN_N=1 \
   CPPAW_CUSOLVER_ACC_CHECK=1 \
   PAWX="../../../bin/nvhpc_cusolver_acc_profile/paw_nvhpc_cusolver_acc_profile.x" \
   make all
   ```
5. If this does not work, the defaults in the parmfile may not be suitable for your system. In this case, copy the default parmfile to the base directory and edit it.
   ```
   cp src/Buildtools/defaultparmfile ./parmfile
   ```
6. After editing the parmfile open `paw_install` and add `-f parmfile` as argument for `paw_build.sh`:
   ```
   src/Buildtools/paw_build.sh -f parmfile -c fast 
   src/Buildtools/paw_build.sh -f parmfile -c fast_parallel
   ```
7. Consult the manual `doc/manual.pdf` for further information.

## License

The CP-PAW code is distributed under the GNU Public License Version 3.
See the LICENSE file.
