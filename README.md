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


# Configuration and Installation instructions


> [!IMPORTANT]
> The installation instructions are for the current version of the CP-PAW code. The instructions may not be applicable to older versions of the code. Installation process is under current development. Please report any issues to the developers.

## Requirements

- fortran compiler (e.g. gfortran, ifort; NVIDIA HPC SDK/nvfortran is supported)
- pkg-config, latexmk, GNU Make (version 4.3 or later)
- bash, cpp, ar
- tex (latex) distribution (e.g. TeX Live) 
- LAPACK, BLAS, FFTW3, MPI (optional), LIBXC (optional)
- optional NVIDIA HPC SDK stack. When `nvfortran` is detected, `./paw_install` also tries `nvhpc_fast` and `nvhpc_fast_parallel` builds using NVPL when available and HPC-X MPI for parallel targets. If CUDA plus NVBLAS are available, optional `nvhpc_nvblas_*` builds link the existing BLAS-3 calls through NVIDIA's GPU BLAS interposition layer. If CUDA plus NVLAMATH are available, optional `nvhpc_nvlamath_*` builds use NVIDIA's LAPACK/cuSOLVER wrapper path. If CUDA plus cuFFTW/cuFFT are available, optional `nvhpc_cufftw_*` builds route CP-PAW's existing FFTW3 calls through cuFFT's FFTW3-compatible wrapper. The experimental `nvhpc_cufft_*` builds require CUDA and use native cuFFT/OpenACC for selected batched 1-D complex FFTs while keeping NVPL FFTW as the fallback. The experimental `nvhpc_cublas_acc_*` builds require CUDA and use OpenACC data regions plus cuBLAS-v2 for selected large complex BLAS-3 kernels. The experimental `nvhpc_cusolver_acc_*` builds require CUDA and use OpenACC data regions plus cuSOLVER-Dn for selected dense real/complex standard and generalized eigensolvers. The combined `nvhpc_gpu_acc_*` builds enable native cuFFT, cuBLAS/OpenACC and cuSOLVER/OpenACC for one-GPU profiling. The opt-in `nvhpc_gpu_all_*` builds link NVPL fallbacks plus cuFFTW, native cuFFT/OpenACC, cuBLAS/OpenACC, cuSOLVER/OpenACC and NVLAMATH in one binary; NVBLAS remains a separate interposition experiment. Set `CPPAW_INSTALL_NVHPC=no` to skip NVIDIA builds, `CPPAW_INSTALL_NVHPC=require` to make them mandatory, `CPPAW_INSTALL_NVBLAS=no` to skip NVBLAS variants, `CPPAW_INSTALL_NVLAMATH=no` to skip NVLAMATH variants, `CPPAW_INSTALL_CUFFTW=no` to skip cuFFTW variants, `CPPAW_INSTALL_CUFFT=no` to skip native cuFFT variants, `CPPAW_INSTALL_GPU_ACC=no` to skip combined GPU variants, `CPPAW_INSTALL_GPU_ALL=yes` to add all-library GPU variants, `CPPAW_INSTALL_CUBLAS_ACC=no` to skip cuBLAS/OpenACC variants, or `CPPAW_INSTALL_CUSOLVER_ACC=no` to skip cuSOLVER/OpenACC variants.
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
   Profiled runs write `cppaw_accel_profile.csv` in serial mode and `cppaw_accel_profile.rankNNNNN.csv` in parallel mode. Set `CPPAW_ACCEL_PROFILE_FILE=<prefix>` to choose another file prefix, or `CPPAW_ACCEL_PROFILE=0` to disable collection at run time. In `nvhpc_cufft_*` and `nvhpc_gpu_acc_*` builds, the native cuFFT path is opt-in: set `CPPAW_CUFFT_ACC=1` to enable it and `CPPAW_CUFFT_ACC_MIN_ELEMENTS=<elements>` to adjust the batched FFT offload threshold; `CPPAW_CUFFT_ACC_3D=1` additionally enables the experimental 3-D cuFFT wrapper for diagnostics. Combined GPU targets use explicit NVHPC GPU memory mode by default (`CPPAW_NVHPC_GPU_MEMORY_MODE=separate`); the `nvhpc_gpu_acc_managed_profile` and `nvhpc_gpu_acc_unified_profile` targets build `-gpu=mem:managed` and `-gpu=mem:unified` variants for device-residency experiments. Set `CPPAW_GPU_RESIDENCY=1` in the residency profile targets to let cuBLAS scalarproducts reuse OpenACC-present operands where projection loops already keep wavefunctions on the device. In `nvhpc_cublas_acc_*` and `nvhpc_gpu_acc_*` builds, `CPPAW_CUBLAS_ACC=0` disables the cuBLAS path at run time and `CPPAW_CUBLAS_ACC_MINFLOP=<flops>` adjusts the offload threshold. In `nvhpc_cusolver_acc_*` and `nvhpc_gpu_acc_*` builds, `CPPAW_CUSOLVER_ACC=0` disables the cuSOLVER path at run time and `CPPAW_CUSOLVER_ACC_MIN_N=<n>` adjusts the dense eigensolver offload threshold; `CPPAW_CUSOLVER_MIN_N=<n>` is accepted as a shorter alias. The default cuBLAS threshold is `1e7`; the Si64 profile favors this over more conservative thresholds. The default cuSOLVER threshold is `256`; smaller Si64 profiling runs can set `CPPAW_CUSOLVER_ACC_MIN_N=1` to force offload. Set `CPPAW_CUSOLVER_ACC_CHECK=1` to validate cuSOLVER eigensolver results by residual and orthonormality before accepting them; failures fall back to the CPU LAPACK path.
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
   For the recommended combined GPU profiling path:
   ```
   cd tests/profile/si64
   PAWX="../../../bin/nvhpc_gpu_acc_profile/paw_nvhpc_gpu_acc_profile.x" \
   make all
   ```
   To scan GPU/NVIDIA-library capabilities and run a short diagnostic matrix:
   ```
   src/Tools/Scripts/paw_gpu_capabilities.sh
   cd tests/profile/si64
   NSTEPS=1 ./run_gpu_exploration.sh
   ```
   To force all native paths for diagnostics, add `CPPAW_CUFFT_ACC=1` and
   `CPPAW_CUSOLVER_ACC_MIN_N=1`.
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
