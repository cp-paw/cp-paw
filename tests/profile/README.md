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

To test the explicit cuBLAS/OpenACC path for large complex `ZGEMM`/`ZHERK`
kernels, build an `nvhpc_cublas_acc_*` target. The default offload threshold is
`CPPAW_CUBLAS_ACC_MINFLOP=1e7`, which includes the projection GEMMs. Set
`CPPAW_CUBLAS_ACC=0` to run the same binary with the CPU/NVPL fallback, or set
`CPPAW_CUBLAS_ACC_MINFLOP=1e8` to keep the projection GEMMs on CPU/NVPL:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cublas_acc_profile_parallel
cd tests/profile/si64
PAWX="mpirun -np 4 ../../../bin/nvhpc_cublas_acc_profile_parallel/ppaw_nvhpc_cublas_acc_profile.x" \
make all
```

The run writes `si64_accel_profile.csv` for serial execution, or
`si64_accel_profile.rankNNNNN.csv` for MPI execution, and prints a compact
summary of the instrumented FFT, BLAS-like and MPI all-to-all regions.

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

Set `RUN_CUFFTW=yes` to add a short cuFFTW comparison suite to the overnight
run.
