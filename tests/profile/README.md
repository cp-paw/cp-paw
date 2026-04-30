# CP-PAW Profiling Cases

These cases are not part of the default test suite. They are intended to
exercise larger kernels for CPU/GPU porting decisions.

Build a profiling executable first, for example:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_profile_parallel
```

Then run the Si64 case with an explicit launcher:

```
cd tests/profile/si64
PAWX="mpirun -np 4 ../../../bin/nvhpc_profile_parallel/ppaw_nvhpc_profile.x" make all
```

To test NVIDIA GPU BLAS interposition for the existing `ZGEMM`, `ZHERK`,
`DGEMM` and `DSYRK` calls, build an `nvhpc_nvblas_*` target, put its script
directory in `PATH` and add `NVBLAS=yes`:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvblas_profile_parallel
cd tests/profile/si64
PATH="../../../bin/nvhpc_nvblas_profile_parallel:${PATH}" \
PAWX="mpirun -np 4 ../../../bin/nvhpc_nvblas_profile_parallel/ppaw_nvhpc_nvblas_profile.x" \
NVBLAS=yes make all
```

To test NVIDIA NVLAMATH for the existing LAPACK calls, build an
`nvhpc_nvlamath_*` target. This enables the NVIDIA HPC SDK LAPACK/cuSOLVER
wrapper path through `-gpu=nvlamath` and keeps the regular NVPL BLAS/FFTW
fallbacks for the rest of the build:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_nvlamath_profile
cd tests/profile/si64
PAWX="../../../bin/nvhpc_nvlamath_profile/paw_nvhpc_nvlamath_profile.x" \
make all
```

To test cuFFTW/cuFFT for the existing FFTW3 calls, build an `nvhpc_cufftw_*`
target. This keeps CP-PAW's current FFT decomposition and uses cuFFT via the
FFTW3-compatible wrapper. It is a correctness and profiling probe, not the
default overnight path:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufftw_profile_parallel
cd tests/profile/si64
PAWX="mpirun -np 4 ../../../bin/nvhpc_cufftw_profile_parallel/ppaw_nvhpc_cufftw_profile.x" \
make all
```

To test the experimental native cuFFT/OpenACC path for batched 1-D complex FFTs,
build an `nvhpc_cufft_*` target. Native cuFFT is opt-in because the Si64
profiling case is dominated by small 1-D FFT batches where copy/setup overheads
outweigh device execution. When enabled, this path uses native cuFFT for
`LIB$FFTC8` when the runtime threshold allows it, but keeps NVPL FFTW as the CPU
fallback for the rest of CP-PAW's FFT work:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cufft_profile_parallel
cd tests/profile/si64
CPPAW_CUFFT_ACC=1 \
PAWX="mpirun -np 4 ../../../bin/nvhpc_cufft_profile_parallel/ppaw_nvhpc_cufft_profile.x" \
make all
```

Set `CPPAW_CUFFT_ACC=1` to enable the native path, and set
`CPPAW_CUFFT_ACC_MIN_ELEMENTS=<elements>` to offload only larger batched FFT
calls.

To profile the combined native GPU paths on one GPU, build an
`nvhpc_gpu_acc_*` target. This enables explicit cuBLAS by default, keeps native
cuFFT opt-in, and uses cuSOLVER only above its default size threshold. The same
binary can selectively force or disable each accelerator path:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_profile
cd tests/profile/si64
NSTEPS=1 RANKS=1 CASES="cpu nvpl nvlamath gpu gpu_force_all gpu_no_cufft gpu_no_cublas gpu_no_cusolver gpu_off" ./run_benchmark.sh
NSTEPS=1 RANKS=8 CASES="cpu nvpl" ./run_benchmark.sh
```

The `cpu` case uses the plain GNU/OpenBLAS/FFTW build (`profile` or
`profile_parallel`) as a pre-HPC-SDK reference. For MPI runs it defaults to the
system `mpirun`; set `CPU_MPIRUN=...` to override it.

To test the explicit cuBLAS/OpenACC path for large complex `ZGEMM`/`ZHERK`
kernels, build an `nvhpc_cublas_acc_*` target. The default offload threshold is
`CPPAW_CUBLAS_ACC_MINFLOP=1e7`, which includes the projection GEMMs and was the
best Si64 threshold in the Spark C86C night run. Set `CPPAW_CUBLAS_ACC=0` to run
the same binary with the CPU/NVPL fallback:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile_parallel
cd tests/profile/si64
PAWX="mpirun -np 4 ../../../bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" \
make all
```

The run writes `si64_accel_profile.csv` for serial execution, or
`si64_accel_profile.rankNNNNN.csv` for MPI execution, and prints a compact
summary of the instrumented FFT, BLAS-like and MPI all-to-all regions.

To test the explicit cuSOLVER/OpenACC dense eigensolver path, build an
`nvhpc_cusolver_acc_*` target. The default offload threshold is
`CPPAW_CUSOLVER_ACC_MIN_N=256`; set it to `1` for the small Si64 profiling
case, or set `CPPAW_CUSOLVER_ACC=0` to run the same binary with the CPU/NVPL
fallback:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_profile
cd tests/profile/si64
NSTEPS=1 RANKS=1 CASES="nvpl cusolver cusolver_off" ./run_benchmark.sh
```

For reproducible comparisons, use the benchmark harness:

```
cd tests/profile/si64
NSTEPS=20 RANKS=4 REPEATS=3 CASES="nvpl cublas cublas_off" ./run_benchmark.sh
```

The harness creates timestamped directories under `tests/profile/si64/runs`,
writes per-run logs and profile CSV files, and emits a `benchmark.tsv` summary
with wall time, instrumented rank-seconds, category timings and final energy.

For a short Nsight Systems trace of the cuBLAS/OpenACC profile binary:

```
cd tests/profile/si64
NSTEPS=1 RANKS=4 ./run_nsys.sh
```

By default the Nsight harness wraps `mpirun` and writes one combined
`nsys_mpi.nsys-rep` report, which is robust for CP-PAW's current
`MPI_ABORT(0)` shutdown path. For runs that finalize MPI normally,
`NSYS_MPI_MODE=per_rank` writes one report per rank.

For an overnight comparison that combines a longer 4-rank run, rank scaling,
a cuBLAS offload-threshold sweep and a short Nsight trace:

```
cd tests/profile/si64
./run_overnight.sh
```

The top-level run directory is written to `runs/latest_overnight`; the combined
benchmark table is `combined_benchmark.tsv`.

Set `RUN_NVLAMATH=yes` to add a short NVLAMATH comparison suite,
`RUN_CUFFTW=yes` to add a short cuFFTW comparison suite, `RUN_CUFFT=yes` to add
a short native cuFFT comparison, `RUN_GPU_ACC=yes` to compare one-rank GPU
against eight-rank CPU/OpenBLAS and CPU/NVPL references, or
`RUN_CUSOLVER=yes` to add the same one-rank cuSOLVER versus eight-rank CPU/NVPL
resource comparison.
