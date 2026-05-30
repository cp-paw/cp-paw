# NVHPC Spark Si64 Benchmark Summary

This note summarizes the Spark C86C Si64 band benchmarks used to choose the
current NVIDIA HPC SDK development defaults. The case is the periodic Si64
profile input with `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`, one GPU
rank for GPU cases, and one-rank plus eight-rank CPU references.

The intent is decision support, not a universal performance claim. Si64 is a
good smoke and orthogonalization/projection case, but larger systems still need
dedicated follow-up runs before promoting any path to production default.

## Run Overview

| Run directory | Scope | Best GPU case | 1-rank CPU references | 8-rank CPU references | Conclusion |
| --- | ---: | ---: | ---: | ---: | --- |
| `si64_bands-nvhpc-standard-20260530-121140` | 6 GPU + CPU | `gpu_resident` 42.86 s | `cpu` 73.35 s, `nvhpc_cpu` 69.73 s | `cpu` 167.41 s, `nvhpc_cpu` 166.76 s | First clear signal that residency helps. |
| `si64_bands-nvhpc-standard-20260530-124102` | CPU-only partial | - | `cpu` 73.32 s, `nvhpc_cpu` 69.45 s | - | Partial run after an invalid case list; CPU context only. |
| `si64_bands-nvhpc-standard-20260530-124708` | Full matrix | `gpu_matmul_conservative` 42.99 s | `cpu` 84.34 s, `nvhpc_cpu` 86.48 s | `cpu` 177.86 s, `nvhpc_cpu` 166.96 s | First full matrix; GPU wins, all-library path is slow. |
| `si64_bands-nvhpc-standard-20260530-135234` | Full matrix, latest | `gpu_resident` 44.64 s | `cpu` 72.91 s, `nvhpc_cpu` 79.66 s | `cpu` 174.13 s, `nvhpc_cpu` 167.66 s | Recommended default direction: residency + explicit cuBLAS. |

The latest run lives at:

```
/home/kuehne88/cp-paw-nvhpc-gpufull/tests/profile/si64/runs/si64_bands-nvhpc-standard-20260530-135234
```

## Full Matrix Comparison

| Case | Previous full matrix | Latest full matrix | Change | Interpretation |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident` | 44.99 s | 44.64 s | -0.8% | Best current default candidate. |
| `gpu_conservative` | 45.31 s | 45.20 s | -0.2% | Very close to the best path. |
| `gpu_addproduct_conservative` | 45.11 s | 45.63 s | +1.2% | Still excellent. |
| `gpu_projection_conservative` | 45.12 s | 45.67 s | +1.2% | Still excellent. |
| `gpu_resident_overlap_conservative` | 45.00 s | 45.92 s | +2.0% | Still excellent. |
| `gpu_resident_projection_conservative` | 44.87 s | 45.95 s | +2.4% | Still excellent. |
| `gpu_resident_addproduct_conservative` | 44.86 s | 46.10 s | +2.8% | Still excellent. |
| `gpu_resident_no_cusolver` | 47.97 s | 46.71 s | -2.6% | cuSOLVER is not the main Si64 lever. |
| `cublas` | 47.72 s | 46.89 s | -1.7% | Good, but residency is better. |
| `gpu_resident_matmul_conservative` | 44.94 s | 46.97 s | +4.5% | Still good. |
| `gpu_matmul_conservative` | 42.99 s | 47.16 s | +9.7% | Previous best was likely run variance or a narrow-case artifact. |
| `cublas_nosync` | 49.19 s | 47.92 s | -2.6% | Diagnostic only until longer correctness runs are available. |
| `gpu_no_cufft` | 45.33 s | 48.12 s | +6.2% | Confirms cuFFT is not decisive for this case. |
| `gpu_no_cusolver` | 53.69 s | 51.21 s | -4.6% | cuSOLVER helps less than residency/cuBLAS. |
| `gpu_off` | 74.71 s | 69.28 s | -7.3% | Same binary without native GPU paths; useful fallback reference. |
| `cufft` | 69.60 s | 74.96 s | +7.7% | Native cuFFT alone is not attractive here. |
| `nvlamath` | 72.21 s | 75.82 s | +5.0% | No Si64 gain. |
| `gpu_no_cublas` | 77.12 s | 79.85 s | +3.5% | cuBLAS is central to the speedup. |
| `nvblas` | 100.10 s | 91.03 s | -9.1% | Still too slow as a default path. |
| `gpu_all` | 161.55 s | 157.59 s | -2.5% | All-library diagnostic path is slow. |
| `cufftw` | 169.50 s | 166.01 s | -2.1% | cuFFTW is not attractive for this case. |
| `gpu_all_3dfft` | 228.71 s | 169.73 s | -25.8% | Improved, but still slow. |
| `gpu_all_off` | 174.37 s | 171.40 s | -1.7% | Slow; useful only as all-library fallback diagnostic. |
| `gpu_all_invbatch_off` | 162.34 s | 179.33 s | +10.5% | Slow; keep as diagnostic only. |

## Current Conclusions

1. Use the residency profile path as the recommended NVHPC GPU profiling path:
   `nvhpc_gpu_acc_residency_profile` and
   `nvhpc_gpu_acc_residency_profile_parallel`.

2. The main win is device residency around wavefunction-heavy regions plus
   explicit cuBLAS. In the latest run, `gpu_resident` is 44.64 s versus
   72.91 s for the one-rank plain CPU reference.

3. Do not make all optional NVIDIA libraries active by default. The
   `gpu_all*`, `cufftw`, `nvblas`, and `nvlamath` cases are valuable diagnostics
   but are slower for this workload.

4. Keep native cuFFT and cuSOLVER threshold-gated. The Si64 result does not
   justify aggressive defaults for either one, although larger generalized
   eigensolver cases may change the cuSOLVER decision.

5. The next implementation target should follow the wavefunction-residency
   question: make data movement around `PSI0`, `PSIM`, `OPSI`, `PROPSI`,
   projection, overlap, and addproduct visible first, then extend the resident
   region only where the profile shows real host/device traffic.

## Recommended Next Benchmark

Use the focused default comparison for routine checks:

```
cd tests/profile/si64
TEST=si64_bands EMPTY_BANDS=1024 NSTEPS=1 RUN_GPU_ALL=no ./run_nvhpc_standard.sh
```

Use the full diagnostic sweep only when comparing library combinations:

```
cd tests/profile/si64
TEST=si64_bands EMPTY_BANDS=1024 NSTEPS=1 RUN_GPU_ALL=yes \
  GPU_CASES="gpu_off gpu_all gpu_all_off gpu_resident gpu_resident_no_cusolver cublas cusolver cufft cufftw nvlamath nvblas gpu_no_cufft gpu_no_cublas gpu_no_cusolver gpu_managed gpu_unified" \
  ./run_nvhpc_standard.sh
```
