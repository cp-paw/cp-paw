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
- optional NVIDIA HPC SDK stack. When `nvfortran` is detected, `./paw_install` also tries `nvhpc_fast` and `nvhpc_fast_parallel` builds using NVPL when available and HPC-X MPI for parallel targets. If CUDA plus NVBLAS are available, optional `nvhpc_nvblas_*` builds link the existing BLAS-3 calls through NVIDIA's GPU BLAS interposition layer. If CUDA plus cuFFTW/cuFFT are available, optional `nvhpc_cufftw_*` builds route CP-PAW's existing FFTW3 calls through cuFFT. The experimental `nvhpc_cublas_acc_*` builds require CUDA and use OpenACC data regions plus cuBLAS-v2 for selected large complex BLAS-3 kernels. Set `CPPAW_INSTALL_NVHPC=no` to skip NVIDIA builds, `CPPAW_INSTALL_NVHPC=require` to make them mandatory, `CPPAW_INSTALL_NVBLAS=no` to skip NVBLAS variants, `CPPAW_INSTALL_CUFFTW=no` to skip cuFFTW variants, or `CPPAW_INSTALL_CUBLAS_ACC=no` to skip cuBLAS/OpenACC variants.
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
   to make those builds mandatory. CUDA-dependent NVBLAS, cuFFTW and cuBLAS/OpenACC variants are only attempted when CUDA is detected, unless requested explicitly with `CPPAW_INSTALL_NVBLAS=require`, `CPPAW_INSTALL_CUFFTW=require` or `CPPAW_INSTALL_CUBLAS_ACC=require`.
   The NVIDIA builds can also be selected directly:
   ```
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_fast_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_fast
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_fast_parallel
   ```
   Profiling builds add CP-PAW hotspot instrumentation for FFT, BLAS/LAPACK-style kernels and MPI transposes:
   ```
   src/Buildtools/paw_build.sh -c profile
   src/Buildtools/paw_build.sh -c profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_profile_parallel
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile
   CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile_parallel
   ```
   Profiled runs write `cppaw_accel_profile.csv` in serial mode and `cppaw_accel_profile.rankNNNNN.csv` in parallel mode. Set `CPPAW_ACCEL_PROFILE_FILE=<prefix>` to choose another file prefix, or `CPPAW_ACCEL_PROFILE=0` to disable collection at run time. In `nvhpc_cublas_acc_*` builds, `CPPAW_CUBLAS_ACC=0` disables the cuBLAS path at run time and `CPPAW_CUBLAS_ACC_MINFLOP=<flops>` adjusts the offload threshold. The default threshold is `1e7`; set `CPPAW_CUBLAS_ACC_MINFLOP=1e8` to keep the long-skinny projection kernels on CPU/NVPL.
   To include profile targets in the installer run:
   ```
   CPPAW_INSTALL_PROFILE=yes ./paw_install
   ```
   The installer compiles with `CPPAW_INSTALL_JOBS=16` by default; set another
   value if your build host needs a smaller or larger parallel make.
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
   To try cuFFTW/cuFFT for CP-PAW's existing FFTW3 calls:
   ```
   cd tests/profile/si64
   PAWX="mpirun -np 4 ../../../bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" \
   make all
   ```
   For explicit cuBLAS/OpenACC profiling without NVBLAS interposition:
   ```
   cd tests/profile/si64
   PAWX="mpirun -np 4 ../../../bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" \
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
