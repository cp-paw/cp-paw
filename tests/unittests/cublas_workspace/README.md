# cuBLAS persistent workspace verification

Build from an existing combined GPU build directory using

```sh
make -f Makefile -f ../../tests/unittests/cublas_workspace/driver.mk cublas-workspace-build
export CPPAW_CUBLAS_ACC=1 CPPAW_CUBLAS_ACC_MINFLOP=1
export CPPAW_CUBLAS_ACC_MATMUL_MINFLOP=1 CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP=1
CPPAW_CUBLAS_FP64_EMULATION=0 ./unit-tests/cublas_workspace.x
CPPAW_CUBLAS_FP64_EMULATION=1 CPPAW_CUBLAS_FP64_TELEMETRY=1 \
  CPPAW_CUBLAS_FP64_STRATEGY=eager CPPAW_CUBLAS_FP64_WORKSPACE_MB=256 \
  ./unit-tests/cublas_workspace.x 1024 used
```

Three alternating DGEMM and ZGEMM calls check CPU reference products, workspace
and telemetry pointer persistence, and freshly written mantissa counters. A
counter of -1 denotes native FP64 fallback and a nonnegative counter denotes
emulation. An unchanged reset value of -2 is unknown, not proof of fallback.
The optional second argument `used` or `fallback` requires that classification
on every call. Emulation without telemetry and a zero user-workspace limit can be
tested by overriding the corresponding environment variables. Matrix sizes are
correctness cases, not timing benchmarks. Builds without the emulation API skip
this test explicitly.

Run the executable under CUDA Compute Sanitizer `--tool memcheck
--error-exitcode 77` to check allocations, copies, and kernel registration. The
Ozaki buffers use the CUDA Runtime C API. OpenACC objects do not require the
additional CUDA Fortran `-cuda` compilation flag, which causes duplicate kernel
registration with NVHPC 26.5 on the tested CUDA 13.2 configuration. The link step
still selects the CUDA runtime. Explicit user requests for CUDA Fortran remain
available through the existing build option.
