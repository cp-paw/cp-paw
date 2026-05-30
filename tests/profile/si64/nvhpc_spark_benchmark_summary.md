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
| `si64_bands-nstep1-1ranks-20260530-211748` | PRO cache smoke | `gpu_resident` 44.08 s | - | - | Cached resident `PRO` beats the host-PRO path for this focused check. |
| `pro-cache-sweep-20260530-213304` | 512/1024/2048 band sweep | `gpu_resident` by a small margin at 1024/3 | Included | Included | GPU residency dominates; full `PRO` cache saves traffic but is near-neutral in wall time. |
| `addpro-cache-split-20260530-231839` | ADDPRO-cache split | `gpu_resident_addpro_host` 45.42 s at 1024/1 | - | 4-rank smoke OK | Adds a diagnostic split between projection cache and `WAVES_ADDPRO` cache reuse. |

The latest full-matrix run lives at:

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

5. The first projector follow-up is now in place: resident `PRO` is built once,
   cached on the GPU and reused by projection plus eligible addproduct calls.
   Keep that full path as the residency default for now because it reduces copy
   estimates substantially, but keep `gpu_resident_addpro_host` in standard
   sweeps because its wall time can be marginally better at 1024-band size.

6. The next implementation target should follow the wavefunction-residency
   question beyond the current orthogonalization envelope: reduce the remaining
   `PSI`/`PROPSI` host/device traffic around projection, overlap and addproduct,
   then retest on larger band/system cases where cache memory and reuse matter
   more than the Si64 smoke.

## Present-Check Smoke

After adding semantic OpenACC present checks, the patched
`nvhpc_gpu_acc_residency_profile` target was rebuilt on Spark C86C and run with
`TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`, one GPU rank:

```
runs/si64_bands-nstep1-1ranks-20260530-192928/gpu_resident/rep01
```

The run completed normally in 45.19 s with final constant energy
302.280854 Ha. The new residency rows reported:

| Profile row | Calls | Estimated copy GB | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_PRESENT_PROJ_PSI` | 2 | 0.0000 | Projection calls inside the resident orthogonalization region reuse `PSI` on device. |
| `ACC_COPY_PROJ_PSI_IN` | 4 | 0.4841 | Projection calls outside that region still copy wavefunctions to the device. |
| `ACC_COPY_ORTHO_PSIM_IN` | 1 | 0.1210 | `PSIM` enters the orthogonalization resident region once. |
| `ACC_COPY_ORTHO_OPSI_IN` | 1 | 0.1210 | `OPSI` enters the orthogonalization resident region once. |
| `ACC_COPY_PROJ_PROPSI_OUT` | 6 | 0.0460 | Projection results are still returned to host for MPI combine and one-center overlap work. |
| `ACC_COPY_ADDPRO_PSI_IO` | 2 | 0.2421 | `WAVES_ADDPRO` still copies its updated wavefunction array in/out. |
| `ACC_COPY_PROJ_PRO_IN` | 384 | 1.0490 | Projector blocks are generated on host and copied for each atom/projection GEMM. |

This supports Peter's concern in a more concrete way: keeping `PSIM`/`OPSI`
resident already works locally in the orthogonalization envelope, but broader
wavefunction residency would need to cover projection calls outside that
envelope, and a separate larger target is moving or caching projector expansion
data (`PRO`) rather than only toggling FFT/LAPACK libraries.

## GPU Projector Expansion Smoke

The follow-up patch moves resident-path projector expansion into OpenACC and
adds a narrow `WAVES_GRAMSCHMIDT` projection/overlap resident region. Spark C86C
builds succeeded for both:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel
```

Smoke run:

```
runs/si64_bands-nstep1-1ranks-20260530-202630
```

| Case | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident` | 54.02 s | 15.7796 GB | 302.280854 Ha | New default path: GPU projector expansion. |
| `gpu_resident_pro_host` | 53.30 s | 16.7391 GB | 302.280854 Ha | Same binary with `CPPAW_GPU_PRO_EXPANSION=0`. |

Key projector rows:

| Profile row | `gpu_resident` | `gpu_resident_pro_host` | Meaning |
| --- | ---: | ---: | --- |
| `ACC_PRESENT_PROJ_PRO` | 384 calls, 0 GB | - | `PRO` is generated and consumed on device. |
| `ACC_COPY_PROJ_PRO_IN` | - | 384 calls, 1.0490 GB | Old host-expanded projector transfer. |
| `ACC_COPY_PROJ_BAREPRO_IN` | 6 calls, 0.0032 GB | - | Bare radial projectors copied once per projection call. |
| `ACC_COPY_PROJ_YLM_IN` | 6 calls, 0.0057 GB | - | Spherical harmonics copied once per projection call. |
| `ACC_COPY_PROJ_EIGR_UPDATE` | 384 calls, 0.0807 GB | - | Structure factor still updated per atom. |
| `ACC_PRESENT_PROJ_PSI` | 4 calls | 4 calls | Broader residency feeds more projection calls with present `PSI`. |
| `ACC_COPY_PROJ_PSI_IN` | 2 calls, 0.2421 GB | 2 calls, 0.2421 GB | Remaining projection calls outside resident wavefunction regions. |
| `ACC_COPY_GRAM_PSI_IN` | 2 calls, 0.2421 GB | 2 calls, 0.2421 GB | New Gramschmidt resident projection/overlap region. |

Conclusion: correctness is preserved and the explicit `PRO` copy drops by about
0.96 GB net in this Si64 smoke. The wall time does not improve yet because the
current OpenACC expansion kernel and per-atom `EIGR` update cost slightly more
than the removed transfer on this small case. Keep GPU projector expansion as a
diagnostic/default in the residency profile for larger-system testing, but the
next optimization should fuse or cache more of the structure-factor/projector
work before claiming a speedup.

## Resident PRO Cache Smoke

The next patch caches the full resident `PRO(NGL,LMNXX,NAT)` block on the GPU,
stores per-atom structure factors in the same cache, and reuses it from both
`WAVES_PROJECTIONS` and eligible `WAVES_ADDPRO` calls. Spark C86C builds
succeeded for both serial and parallel residency targets:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel
```

Focused one-rank smoke:

```
runs/si64_bands-nstep1-1ranks-20260530-211748
```

| Case | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident` | 44.08 s | 15.3405 GB | 302.280854 Ha | Cached GPU `PRO` plus cached addproduct path. |
| `gpu_resident_pro_host` | 47.24 s | 16.4970 GB | 302.280854 Ha | Same binary with `CPPAW_GPU_PRO_EXPANSION=0`. |

Key cache rows from the one-rank smoke:

| Profile row | `gpu_resident` | `gpu_resident_pro_host` | Meaning |
| --- | ---: | ---: | --- |
| `ACC_BUILD_PRO_CACHE` | 1 call, 0.0491 s, 0.1897 GB | - | Full resident projector cache is built once for this geometry/grid. |
| `ACC_PRESENT_PRO_CACHE_REUSE` | 7 calls | - | Later projection/addproduct uses reuse the cache. |
| `ACC_PRESENT_PROJ_PRO_CACHE` | 6 calls, 0 GB | - | Projection consumes cached `PRO` on device. |
| `ACC_COPY_PROJ_PRO_IN` | - | 384 calls, 1.0490 GB | Old host-expanded projector transfer. |
| `ACC_PRESENT_ADDPRO_PRO_CACHE` | 2 calls, 0 GB | - | `WAVES_ADDPRO` consumes cached `PRO` on device. |
| `CUBLAS_ZGEMM_ADDPRO_CACHE` | 128 calls, 0.3719 s | - | Addproduct update is routed through the present-cache cuBLAS path. |

Parallel 4-rank smoke:

```
runs/si64_bands-nstep1-4ranks-20260530-212047
```

| Case | Ranks | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident` | 4 | 103.72 s | 26.1007 GB | 302.280854 Ha | Parallel binary exercises the same cache path; performance is not a target because four ranks share one GPU. |

Each MPI rank built the cache once and reused it seven times; the per-rank
`ACC_BUILD_PRO_CACHE` cost was 0.063-0.070 s. This makes the cache path a useful
baseline for larger systems and for the next residency step, but not yet proof
that the whole PAW wavefunction path should stay permanently resident.

## Larger PRO Cache Sweep

The focused sweep below used the same Spark C86C branch before adding the
separate `ADDPRO` cache switch:

```
runs/pro-cache-sweep-20260530-213304
```

Completed cases were `EMPTY_BANDS=512,1024` with `NSTEPS=1,3`, plus
`EMPTY_BANDS=2048` with `NSTEPS=1`. The planned `2048,NSTEPS=3` point was
stopped after the 2048/1 CPU references because the run time was no longer
reasonable for an interactive pass; keep it as a separate night-run candidate.

| Empty bands | NSTEPS | `gpu_resident` | `gpu_resident_pro_host` | `gpu_off` | 1-rank CPU best | 8-rank CPU best | Interpretation |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 512 | 1 | 12.19 s | 11.07 s | 19.30 s | 19.37 s | 32.37 s | GPU residency helps; full cache not yet a wall-time win. |
| 512 | 3 | 20.80 s | 20.79 s | 38.14 s | 38.61 s | 49.82 s | Cache and host-PRO are tied; residency dominates. |
| 1024 | 1 | 46.35 s | 45.98 s | 75.07 s | 69.64 s | 167.84 s | GPU residency clearly wins; `PRO` cache is near-neutral. |
| 1024 | 3 | 69.67 s | 70.08 s | 134.36 s | 136.53 s | 245.48 s | Full cache is slightly faster and saves about 2.25 GB estimated copies. |
| 2048 | 1 | 281.70 s | 281.02 s | 393.06 s | 386.81 s | 1130.65 s | GPU residency remains the main win; `PRO` cache is wall-time neutral. |

The 2048/1 profile makes the next bottleneck clearer: `gpu_resident` records
only about 24.4 instrumented rank-seconds inside 281.7 wall seconds. BLAS,
LAPACK and FFT offload explain the CPU/GPU difference, but much of the remaining
wall time is still outside the current accelerator timers. That points the next
implementation pass toward broader PAW/wavefunction-region instrumentation and
residency, not toward adding another NVIDIA library first.

## ADDPRO Cache Split

The follow-up patch adds `CPPAW_GPU_ADDPRO_CACHE=0` and the benchmark case
`gpu_resident_addpro_host`. This keeps the GPU `PRO` cache for projections but
routes `WAVES_ADDPRO` through the previous addproduct path, so projection-cache
effects and addproduct-cache effects can be compared separately.

```
runs/addpro-cache-split-20260530-231839
```

| Empty bands | Ranks | `gpu_resident` | `gpu_resident_addpro_host` | `gpu_resident_pro_host` | Energy |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 1 | 11.19 s | 11.99 s | 11.14 s | 302.280854 Ha |
| 1024 | 1 | 45.90 s | 45.42 s | 45.56 s | 302.280854 Ha |
| 512 | 4 | - | 22.17 s | - | 302.280854 Ha |

The split case is correct and useful, but it does not justify changing the
default yet. At 1024 bands it is slightly faster, while the full cache path
keeps much lower estimated `WAVES_ADDPRO` transfer volume
(`ACC_COPY_ADDPRO_PROPSI_IN` plus `CUBLAS_ZGEMM_ADDPRO_CACHE` instead of the
large generic `CUBLAS_ZGEMM_ADDPRODUCT` copy estimate). Keep the full cache as
the default and include `gpu_resident_addpro_host` in future standard/large-band
sweeps.

## PAW Envelope Profiling

The next profiling patch adds coarse `PAW_*` envelope rows around the
wavefunction and PAW regions that were previously mostly invisible in the
accelerator CSV: `WAVES$HPSI`, its `VPSI`/`HPROJ`/`ADDPRO` sections,
`WAVES_ADDPRO`, `WAVES_OPSI`, `WAVES_OPROJ`, `WAVES_ADDOPSI`,
`WAVES_1COVERLAP`, `WAVES$RHO`, `WAVES_DENSITY`, and `WAVES$SPHERE`.
These rows intentionally double-count lower-level BLAS/FFT rows; their purpose
is to show where the remaining wall time sits before moving more data and
operations into resident GPU regions.

A Spark smoke test with `EMPTY_BANDS=512,NSTEPS=1` verifies the new rows.
`coverage` remains the non-envelope coverage; `paw_s` is nested diagnostic
time and can overlap with BLAS/FFT/LAPACK rows.

| Run | Ranks | Wall time | Coverage | `paw_s` | Energy | Dominant new envelope rows |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `paw-envelope-smoke-20260530-234323-1rank` | 1 | 11.32 s | 40.33% | 5.13 s | 302.280854 Ha | `PAW_1COVERLAP_TOTAL`, `PAW_HPSI_TOTAL`, `PAW_HPSI_VPSI`, `PAW_RHO_TOTAL` |
| `paw-envelope-smoke-20260530-234349-4rank` | 4 | 22.15 s | 19.72% | 13.16 s | 302.280854 Ha | `PAW_HPSI_TOTAL`, `PAW_1COVERLAP_TOTAL`, `PAW_HPSI_VPSI`, `PAW_ADDPRO_TOTAL` |

## One-Center Overlap Offload

The next residency patch moves the dense `WAVES_1COVERLAP` contraction to
cuBLAS when `CPPVAR_CUBLAS_ACC` is available. This is intentionally separate
from the existing wavefunction-overlap residency path because the one-center
projection arrays are packed from the PAW projector layout first. The first
Spark smoke result was correct but slower (`11.89 s` for opt-in offload versus
`10.96 s` with `CPPAW_GPU_1COVERLAP=0` at `EMPTY_BANDS=512,NSTEPS=1`), and a
larger `2048/1` probe showed the GPU idle while the process spun on CPU. Keep
this path opt-in for now. Use `CPPAW_GPU_1COVERLAP=1` or
`CPPAW_CUBLAS_ACC_1COVERLAP=1` to enable it explicitly; use the benchmark case
`gpu_resident_1coverlap`. The explicit host diagnostic remains
`gpu_resident_1coverlap_host`.

After switching the default off again, the safe residency path rebuilds and runs
normally. The 512/1 smoke keeps the same energy for `gpu_resident` and
`gpu_resident_1coverlap_host` (`302.280854 Ha`). The larger default-off 2048/1
check completed as
`onecenter-defaultoff-2048-nstep1-20260531-012732`: `288.08 s` wall time,
`302.280854 Ha`, and `PAW_1COVERLAP_TOTAL=13.57 s`. No `CUBLAS_ZGEMM_1COV_*`
rows are emitted unless the opt-in keyword is set. The explicit opt-in smoke
`onecenter-optin-smoke-20260531-013314` completed correctly at 512/1
(`12.40 s`, `302.280854 Ha`) and emitted the expected `CUBLAS_ZGEMM_1COV_*`
profile rows.

For the expensive unresolved point, use the dedicated night-run wrapper. It
defaults to the GPU cases only so the run is not dominated by the known slow
8-rank CPU reference; add `RUN_CPU_REFERENCES=yes` when CPU reference numbers
are explicitly needed.

```
cd tests/profile/si64
EMPTY_BANDS=2048 NSTEPS=3 ./run_gap_profile_night.sh
```

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
