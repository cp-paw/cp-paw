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
> The NVIDIA accelerator paths on this branch are
> experimental. CPU fallbacks remain available, and optional accelerator
> features are either selected by a dedicated build target or guarded by a
> run-time switch. Validate energies, forces, stress, and parallel parity for
> the intended system before using these paths for production calculations.

This branch remains fully usable without an NVIDIA GPU, CUDA, or the NVIDIA
HPC SDK. The standard `dbg`, `fast`, and `fast_parallel` binaries retain the
main-branch dependency set and are built first. Use
`CPPAW_INSTALL_NVHPC=no ./paw_install` for a strictly GPU-free installation.
CI builds serial and MPI GNU configurations without NVIDIA libraries and checks their dynamic dependencies.

The build system detects the NVIDIA HPC SDK, CUDA, and individual libraries
instead of assuming that every NVIDIA installation provides the same stack.
With the installer's default `auto` policy, unavailable optional targets are
skipped. Set an installer option to `require` to turn a missing dependency or
failed build into an installation error.

### NVIDIA build profiles

The normal user interface has three build profiles. Library-by-library targets
remain available for development and reproducible profiling, but are not
separate recommended CP-PAW variants.

| Profile | Purpose | Main libraries |
| --- | --- | --- |
| `nvhpc_fast*` | NVIDIA CPU build | NVPL BLAS, LAPACK, and FFTW when available |
| `nvhpc_gpu_fast*` | Recommended NVIDIA GPU build | cuBLAS, cuSOLVER, native cuFFT, OpenACC residency |
| `nvhpc_gpu_profile*` | Instrumented form of the recommended GPU build | same as `nvhpc_gpu_fast*` plus accelerator telemetry |

The GPU executable provides three run-time modes through `CPPAW_GPU_MODE`:

- `resident` is the default and keeps wavefunctions, projections, and selected
  work arrays on the GPU across neighboring operations. On one MPI rank it
  also uses the native full-grid 3-D cuFFT path with device-side mapping.
- `transfer` keeps cuBLAS and cuSOLVER available but disables the broad
  residency and native 3-D FFT defaults. It is the main comparison and
  compatibility mode.
- `off` disables the explicit cuBLAS, cuSOLVER, and native cuFFT paths in the
  same executable. Use `nvhpc_fast*` for a clean CPU performance reference.

```sh
# CPU and recommended GPU release builds
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast -j16 -z
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_fast -j16 -z
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_fast_parallel -j16 -z

# Same GPU profile with instrumentation
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_profile -j16 -z
```

The `resident` defaults use a cuBLAS offload threshold of `1e7`, cuSOLVER for
supported problems of size 256 or larger, and Gram-Cholesky at size 4096 or
larger. CUDA 13 FP64 Ozaki emulation remains disabled until
`CPPAW_CUBLAS_FP64_EMULATION=1` is set.

The profile names describe user intent rather than every linked library
combination. Use the profiling harness to compare them for the intended
system and hardware.

Advanced targets include `nvhpc_gpu_acc_*` for the unconfigured combined
libraries, `nvhpc_gpu_acc_residency_profile*` for the earlier focused-residency
experiments, and individual cuBLAS, cuSOLVER, cuFFT, cuFFTW, NVLAMATH, and
NVBLAS targets. `nvhpc_gpu_all_*` is an integration test that additionally
links cuFFTW and NVLAMATH. It is an integration configuration, not a
fastest-library preset.
NVBLAS remains separate because it interposes the host BLAS interface.

Skala development is separate in the [`skala`](https://github.com/cp-paw/cp-paw/tree/skala) branch, based on `main`.
The preserved [`nvhpc-skala`](https://github.com/cp-paw/cp-paw/tree/nvhpc-skala) branch retains the complete combined implementation, including Skala GPU kernels and their residency interfaces.

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

### Profiling and validation

Profiling builds record FFT, BLAS/LAPACK, accelerator, MPI transpose, transfer,
and residency data in `cppaw_accel_profile*.csv`. The Si64 harness
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
   This always builds the conventional `dbg`, `fast`, and `fast_parallel`
   binaries. To disable discovery and construction of every NVIDIA-specific
   variant explicitly, use:
   ```
   CPPAW_INSTALL_NVHPC=no ./paw_install
   ```
   On systems where the NVIDIA HPC SDK is available, the installer will additionally try CPU/NVPL builds. Use:
   ```
   CPPAW_INSTALL_NVHPC=require ./paw_install
   ```
   to make those builds mandatory. The installer adds `nvhpc_gpu_fast*` only
   when CUDA, cuFFT, cuBLAS, and cuSOLVER are detected. Set
   `CPPAW_INSTALL_GPU_ACC=require` to make that profile mandatory. Individual
   and all-library targets are opt-in diagnostics; their installer switches
   are described in `tests/profile/README.md`.
   With `CPPAW_INSTALL_PROFILE=yes`, the installer also builds
   `nvhpc_gpu_profile*`. The legacy focused-residency profile is skipped by
   default and remains available with
   `CPPAW_INSTALL_GPU_RESIDENCY_PROFILE=yes`.
   The NVIDIA builds can also be selected directly:
   ```
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_fast_parallel
   ```
   Profiling builds add CP-PAW hotspot instrumentation for FFT, BLAS/LAPACK-style kernels and MPI transposes:
   ```
   src/Buildtools/paw_build.sh -c profile
   src/Buildtools/paw_build.sh -c profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_profile_parallel
   ```
   Profiled runs write `cppaw_accel_profile.csv` in serial mode and
   `cppaw_accel_profile.rankNNNNN.csv` in parallel mode. Set
   `CPPAW_ACCEL_PROFILE_FILE=<prefix>` to choose another prefix, or
   `CPPAW_ACCEL_PROFILE=0` to disable collection. Fine-grained accelerator
   switches, thresholds, memory-mode experiments, and diagnostic build targets
   are documented in `tests/profile/README.md`.
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
   PAWX="../../../bin/nvhpc_gpu_profile/paw_nvhpc_gpu_profile.x" \
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
   cuFFT, cuBLAS, cuSOLVER, and a usable host numerical stack it points routine
   checks at `gpu_recommended`, `gpu_transfer`, and `gpu_off`. The
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
