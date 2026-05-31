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
| `si64_bands-nvhpc-standard-20260531-123955` | Focused standard refresh | `gpu_resident_addpro_host` 43.06 s | `cpu` 77.85 s, `nvhpc_cpu` 75.61 s | `cpu` 167.40 s, `nvhpc_cpu` 167.70 s | Confirms the residency path remains the useful GPU direction on Spark. |
| `addpro-profile-contexts512-20260531-172542` | ADDPRO context profiling | `gpu_resident` 6.66 s at 512/1 | - | `gpu_resident` 9.52 s at 512/4 | Splits ADDPRO copies into `HPSI` and `OPSI` and corrects the `PSI` row to input/output accounting. |
| `fresh-energy-guard-sweep-20260531-174614` | Energy-guard refresh | `gpu_resident_orthox` 13.39 s at 1024/1 | `cpu` 70.83 s, `nvhpc_cpu` 71.85 s | `cpu` 168.36 s, `nvhpc_cpu` 167.20 s | Confirms `WAVES_ORTHO_X` workspace residency is now the best default inside the residency profile. |
| `gram-profile-contexts-20260531-180922` | Gram context profiling | `gpu_resident` 6.94 s at 512/1 | - | `gpu_resident` 9.15 s at 512/4 | Splits the initial Gram wavefunction copy into `PSI0` and `PSIM` rows. |
| `opsi-build-residency-20260531-164302` | OPSI build-residency diagnostic | `gpu_resident`/`gpu_resident_opsi` tied for Si64 | - | `gpu_resident`/`gpu_resident_opsi` tied at 512/4 | Adds an opt-in non-superwave OPSI residency switch; Si64 is a superwave case, so the guard correctly leaves it unchanged. |
| `superwave-opsi-hostscale-20260531-*` | Superwave overlap residency | `gpu_resident`/`gpu_resident_opsi` tied and energy-valid | - | `gpu_resident_opsi` 9.15 s at 512/4 | Makes the superwave inversion overlap term resident; keeps superwave OPSI host-built/host-scaled before entering device residency. |
| `addpro-context-cache-20260531-*` | ADDPRO context-cache controls | `gpu_resident_opsi` 40.03 s at 2048/1 | - | `gpu_resident_addpro_hpsi_host` 9.30 s at 512/4 | Adds independent HPSI/OPSI `WAVES_ADDPRO` cache switches; all cases remain energy-valid, but timings are neutral/noisy. |

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

6. The one-center overlap GPU-pack path removes the previous large
   `WAVES_1COVERLAP` bottleneck. The next large-band hotspot became the initial
   Gram-Schmidt solve inside `WAVES_ORTHO_Y_C`; the Cholesky follow-up below is
   the first direct fix for that path.

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

The first one-center residency patch moved the dense `WAVES_1COVERLAP`
contraction to cuBLAS when `CPPVAR_CUBLAS_ACC` is available. That version was
correct but slower (`11.89 s` for opt-in offload versus `10.96 s` with
`CPPAW_GPU_1COVERLAP=0` at `EMPTY_BANDS=512,NSTEPS=1`), and a larger `2048/1`
probe showed the GPU idle while the process spun on CPU. The follow-up patch
keeps the one-center path but packs the flattened cuBLAS input matrices on the
GPU (`ACC_PACK_CUBLAS_1COV`) before the two `ZGEMM` contractions. With that
fix, the path is enabled by default in residency-profile builds. Use
`CPPAW_GPU_1COVERLAP=0` or `CPPAW_CUBLAS_ACC_1COVERLAP=0` to disable it; the
explicit host diagnostic remains `gpu_resident_1coverlap_host`.

Historical default-off checks showed that the safe residency path rebuilt and ran
normally. The 512/1 smoke keeps the same energy for `gpu_resident` and
`gpu_resident_1coverlap_host` (`302.280854 Ha`). The larger default-off 2048/1
check completed as
`onecenter-defaultoff-2048-nstep1-20260531-012732`: `288.08 s` wall time,
`302.280854 Ha`, and `PAW_1COVERLAP_TOTAL=13.57 s`. No `CUBLAS_ZGEMM_1COV_*`
rows are emitted unless the opt-in keyword is set. The explicit opt-in smoke
`onecenter-optin-smoke-20260531-013314` completed correctly at 512/1
(`12.40 s`, `302.280854 Ha`) and emitted the expected `CUBLAS_ZGEMM_1COV_*`
profile rows.

GPU-pack follow-up smokes:

```
runs/onecoverlap-gpupack-smoke512-20260531-133923
runs/onecoverlap-gpupack-smoke1024-20260531-134008
runs/onecoverlap-gpupack-smoke2048-20260531-134158
runs/onecoverlap-gpupack-final1024-20260531-135518
```

| Case | Empty bands | Wall time | `PAW_1COVERLAP_TOTAL` | Final energy | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| Host contraction | 512 | 11.59 s | 1.0452 s | 302.280854 Ha | Previous default path. |
| GPU-pack cuBLAS contraction | 512 | 10.62 s | 0.0566 s | 302.280854 Ha | New path is faster even in the small smoke. |
| Host contraction | 1024 | 45.54 s | 3.4300 s | 302.280854 Ha | Previous default path. |
| GPU-pack cuBLAS contraction | 1024 | 39.49 s | 0.1472 s | 302.280854 Ha | Best Si64 1024-band signal so far. |
| Host contraction | 2048 | 282.69 s | 12.2053 s | 302.280854 Ha | Previous large-band problem case. |
| GPU-pack cuBLAS contraction | 2048 | 269.85 s | 0.4511 s | 302.280854 Ha | Contract bottleneck is removed; total wall time still has other large-band costs. |

After enabling the GPU-pack path by default, the final 1024-band smoke reported
`gpu_resident` at 41.69 s versus 45.22 s for
`gpu_resident_1coverlap_host`, both with final energy 302.280854 Ha.

For the expensive unresolved point, use the dedicated night-run wrapper. It
defaults to the GPU cases only so the run is not dominated by the known slow
8-rank CPU reference; add `RUN_CPU_REFERENCES=yes` when CPU reference numbers
are explicitly needed.

```
cd tests/profile/si64
EMPTY_BANDS=2048 NSTEPS=3 ./run_gap_profile_night.sh
```

## Inner PAW Split Profiling

The next profiling pass splits the coarse PAW envelope rows into actionable
subregions without changing the numerical path. The new rows distinguish
`PAW_OPSI_OPROJ` from `PAW_OPSI_ADDPRO`, split `WAVES_VPSI` into
`PAW_VPSI_FFT_GTOR`, `PAW_VPSI_POT`, `PAW_VPSI_FFT_RTOG`,
`PAW_VPSI_KIN`, and optional `PAW_VPSI_BUCKET`, and split
`WAVES_DENSITY`/`WAVES$RHO` into FFT, kinetic-density, accumulation, combine,
and spin-conversion rows. `PAW_DENSITY_KIN_*` rows only appear when kinetic
density is requested.

Spark C86C validation after the split:

| Run | Empty bands | Ranks | Wall time | `paw_s` | `blas_s` | `fft_s` | Energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `si64_bands-nstep1-1ranks-20260531-091623` | 512 | 1 | 12.81 s | 7.03 s | 3.19 s | 2.14 s | 302.280854 Ha |
| `si64_bands-nstep1-1ranks-20260531-091705` | 2048 | 1 | 282.52 s | 32.17 s | 16.35 s | 7.04 s | 302.280854 Ha |
| `si64_bands-nstep1-4ranks-20260531-092220` | 512 | 4 | 22.12 s | 17.41 s | 10.20 s | 4.63 s | 302.280854 Ha |

The 2048/1 split confirms the current direction. `PAW_VPSI_TOTAL=2.4375 s`
is almost entirely the two FFT phases
(`PAW_VPSI_FFT_GTOR=1.2045 s`,
`PAW_VPSI_FFT_RTOG=1.1913 s`), while the real-space potential multiply is only
`0.0280 s` and kinetic addition is `0.0137 s`. `PAW_RHO_TOTAL=1.6100 s` is
likewise almost entirely `PAW_DENSITY_FFT=1.5336 s`; the density accumulation is
only `0.0053 s`. `PAW_OPSI_TOTAL=0.3631 s` is dominated by
`PAW_OPSI_ADDPRO=0.3571 s`, and the small one-center `OPROJ` contraction is not
a bottleneck.

Practical conclusion: for this Si64 stress input, accelerating the scalar
real-space loops inside `VPSI` or `RHO` would not move wall time. If we keep more
PAW data on the GPU, the value comes from making FFT and projection phases
consume resident wavefunction/projector data instead of copying around them. The
larger unresolved algorithmic hotspot remains the dense Gram/overlap side:
`PAW_1COVERLAP_TOTAL=12.8553 s` at 2048/1, plus the heavy cuBLAS scalar-product
and matmul rows. That is the better next design target than another local loop
offload.

## Orthogonalization and 1COVERLAP Phase Profiling

The next profiling patch splits the orthogonalization step into context rows and
breaks `WAVES_1COVERLAP` into pack, contraction, and superwave-unravel phases.
This is a diagnostic-only change; it does not alter the numerical path or any
NVHPC defaults.

Spark C86C validation after the split:

| Run | Empty bands | Ranks | Wall time | Energy | Main new signal |
| --- | ---: | ---: | ---: | ---: | --- |
| `si64_bands-nstep1-1ranks-20260531-120036` | 512 | 1 | 11.09 s | 302.280854 Ha | `PAW_1COVERLAP_CONTRACT=1.0021 s` of `1.0320 s`; `PAW_ORTHO_1COVERLAP=0.6184 s`. |
| `si64_bands-nstep1-4ranks-20260531-120104` | 512 | 4 | 21.93 s | 302.280854 Ha | `PAW_ORTHO_SOLVE=2.3215 s`, `PAW_ORTHO_ADDOPROJ=2.1621 s`, `PAW_1COVERLAP_CONTRACT=1.6548 s` rank-summed. |
| `si64_bands-nstep1-1ranks-20260531-120141` | 2048 | 1 | 280.99 s | 302.280854 Ha | `PAW_ORTHO_TOTAL=21.5684 s`, with `PAW_ORTHO_SOLVE=8.6428 s` and `PAW_ORTHO_1COVERLAP=7.5910 s`. |

The 2048/1 result sharpens the next implementation choice. The one-center
overlap bottleneck is not packing or superwave expansion:
`PAW_1COVERLAP_PACK=0.1396 s`,
`PAW_1COVERLAP_UNRAVEL=0.0149 s`, but
`PAW_1COVERLAP_CONTRACT=12.3631 s` out of
`PAW_1COVERLAP_TOTAL=12.5175 s`. The earlier opt-in cuBLAS 1C path should
therefore not be replaced by more host packing work; it needs either a robust
contract-only GPU implementation or an algorithmic change that avoids repeating
the dense contraction in the current form. Orthogonalization also exposes a
second CPU-side target at large band count: the overlap solve/update phase
(`PAW_ORTHO_SOLVE`) is now comparable to the 1C contribution.

## Standard Refresh and MATMUL Residency Split

The 2026-05-31 Spark refresh reran the focused standard comparison on the
current `cp-paw-nvhpc` baseline with `TEST=si64_bands`, `EMPTY_BANDS=1024`,
`NSTEPS=1`, one GPU rank, and one-/eight-rank CPU references:

```
runs/si64_bands-nvhpc-standard-20260531-123955
```

| Case | Ranks | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident` | 1 | 43.70 s | 15.3405 GB | 302.280854 Ha | Main residency default remains fast and correct. |
| `gpu_resident_addpro_host` | 1 | 43.06 s | 15.4480 GB | 302.280854 Ha | Slightly fastest in this run; still kept as a diagnostic split rather than the default. |
| `gpu_resident_invbatch_off` | 1 | 48.54 s | 432.9138 GB | 302.280854 Ha | Confirms inverse-batch residency should stay enabled. |
| `gpu_resident_no_cusolver` | 1 | 48.94 s | 15.2439 GB | 302.280854 Ha | cuSOLVER is useful but not the main Si64 lever. |
| `cpu` | 1 | 77.85 s | - | 302.280854 Ha | Plain CPU reference. |
| `nvhpc_cpu` | 1 | 75.61 s | - | 302.280854 Ha | NVPL/NVHPC CPU reference. |
| `cpu` | 8 | 167.40 s | - | 302.280854 Ha | Resource comparison; not efficient on the Spark CPU state used here. |
| `nvhpc_cpu` | 8 | 167.70 s | - | 302.280854 Ha | NVPL does not rescue this 8-rank Spark comparison. |

The follow-up patch adds resident-aware data regions and granular copy/present
profiling around the generic `LIB$MATMUL` cuBLAS paths. It also fixes NVHPC root
detection when only `nvfortran` from `$NVHPC_ROOT/compilers/bin` is visible in
`PATH`. Spark builds succeeded for both serial and parallel residency targets
with `NVHPC_ROOT` unset.

Post-patch smoke:

```
runs/matmul-residency-smoke-20260531-125804
```

| Case | Ranks | Wall time | Total copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: |
| `gpu_resident` | 1 | 46.02 s | 15.3405 GB | 302.280854 Ha |
| `gpu_resident_addpro_host` | 1 | 45.56 s | 15.4480 GB | 302.280854 Ha |

The old aggregate `ACC_COPY_CUBLAS_ZGEMM_MAT` row is now split by operand in
the residency path:

| Profile row | Calls | Estimated copy GB | Meaning |
| --- | ---: | ---: | --- |
| `ACC_COPY_ZGEMM_MAT_A_IN` | 64 | 7.7462 | Dominant remaining generic complex MATMUL transfer. |
| `ACC_COPY_ZGEMM_MAT_B_IN` | 64 | 0.0077 | Small right-hand operand. |
| `ACC_COPY_ZGEMM_MAT_C_OUT` | 64 | 0.1748 | Small result relative to the left operand. |
| `ACC_COPY_DGEMM_MAT_A_IN` | 58 | 0.6158 | Real MATMUL left operand. |
| `ACC_COPY_DGEMM_MAT_B_IN` | 58 | 0.6158 | Real MATMUL right operand. |
| `ACC_COPY_DGEMM_MAT_C_OUT` | 58 | 0.6158 | Real MATMUL result. |

Conclusion: the major unresolved copy target in this Si64 profile is no longer
ambiguous. The complex generic MATMUL shape is `13133 x 576 x 13`, and almost
all of its transfer volume is the left operand. The next useful implementation
step is to identify that caller's producer and either keep that operand resident
or route it through a more semantic present-data path instead of treating it as
an opaque `LIB$MATMUL` temporary.

## Force PSI Residency

The dominant complex MATMUL transfer from the previous section was traced to
`WAVES_DEDPRO`:

```
CALL LIB$MATMULC8(NGL,NDIM*NBH,LMNX,PSI,DEDPROJ1,DEDPRO)
```

For the Si64 band profile this is the `13133 x 576 x 13` shape. The follow-up
patch keeps `THIS%PSI0` resident across the per-atom `WAVES$FORCE` loop, so the
64 per-atom `WAVES_DEDPRO` calls can reuse the same left MATMUL operand. The
path is enabled by default in residency-profile builds and can be disabled with
`CPPAW_GPU_FORCE_PSI_RESIDENCY=0`.

Spark C86C builds succeeded for both:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel
```

Smoke run:

```
runs/forcepsi-residency-smoke-20260531-130956
```

| Case | Ranks | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident` | 1 | 42.85 s | 7.7153 GB | 302.280854 Ha | New default path with force-loop `THIS%PSI0` residency. |
| `gpu_resident_forcepsi_host` | 1 | 45.56 s | 15.3405 GB | 302.280854 Ha | Same binary with `CPPAW_GPU_FORCE_PSI_RESIDENCY=0`. |

Key profile rows:

| Profile row | `gpu_resident` | `gpu_resident_forcepsi_host` | Meaning |
| --- | ---: | ---: | --- |
| `ACC_COPY_FORCE_PSI0_IN` | 1 call, 0.1210 GB | - | One copy into the force-loop resident region. |
| `ACC_PRESENT_ZGEMM_MAT_A` | 64 calls | - | `WAVES_DEDPRO` reuses the resident `PSI0` operand. |
| `ACC_COPY_ZGEMM_MAT_A_IN` | - | 64 calls, 7.7462 GB | Previous per-atom copy of the same left operand. |
| `ACC_COPY_ZGEMM_MAT_B_IN` | 64 calls, 0.0077 GB | 64 calls, 0.0077 GB | Small per-atom derivative-projector operand. |
| `ACC_COPY_ZGEMM_MAT_C_OUT` | 64 calls, 0.1748 GB | 64 calls, 0.1748 GB | Per-atom `DEDPRO` output. |

Conclusion: this is the first broad force-side wavefunction-residency win. It
does not solve the remaining DGEMM MATMUL traffic or one-center overlap
contraction, but it removes the largest previously ambiguous complex MATMUL copy
without changing the energy and with a positive one-step wall-time signal on
Spark.

## Orthogonalization Constant Residency Diagnostic

The remaining real generic MATMUL traffic is mainly inside `WAVES_ORTHO_X`.
A prototype kept only the constant `CHICHI` and `U` inputs resident across the
safe-orthogonalization loop, while leaving host-updated temporaries on the
existing copy-back path. This reduced DGEMM input-copy volume but did not improve
Spark wall time, so the path is kept as an opt-in diagnostic instead of a
residency-profile default.

Spark C86C validation:

```
runs/orthoconst-residency-repeat-20260531-132042
runs/orthoconst-residency-bands1536-20260531-132557
runs/orthoconst-final-smoke-20260531-133304
```

| Case | Bands setting | Repeats | Wall time | Total copy estimate | Final energy | Interpretation |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Orthogonalization constants resident | `EMPTY_BANDS=1024` | 3 | 44.50 s avg | 7.1208 GB | 302.280854 Ha | Saves about 0.59 GB copy but is slower than the default in this smoke. |
| Default per-call DGEMM input copies | `EMPTY_BANDS=1024` | 3 | 43.50 s avg | 7.7153 GB | 302.280854 Ha | Better wall time on Spark despite more copy volume. |
| Orthogonalization constants resident | `EMPTY_BANDS=1536` | 1 | 129.54 s | 12.4526 GB | 302.280854 Ha | Saves about 1.44 GB copy in the larger band smoke. |
| Default per-call DGEMM input copies | `EMPTY_BANDS=1536` | 1 | 123.72 s | 13.8925 GB | 302.280854 Ha | Still faster on Spark, so default remains off. |

The diagnostic is exposed as `gpu_resident_orthoconst` and can also be enabled
directly with `CPPAW_GPU_ORTHO_CONST_RESIDENCY=1`. This keeps the implementation
available for larger systems or different interconnect/GPU-memory behavior
without penalizing the current recommended Spark path.

The final one-repeat smoke with the finished harness semantics reported
`gpu_resident` at 42.68 s and 7.7153 GB copy versus `gpu_resident_orthoconst`
at 45.08 s and 7.1208 GB copy, both with final energy 302.280854 Ha.

## ETOT and Gram-Schmidt Profiling

After the one-center GPU-pack path, the 2048-band profile still showed a large
`PHASE_ETOT_WAVES` block with too little internal structure. The follow-up
profiling patch splits `WAVES$ETOT` into `PAW_ETOT_*` rows and further splits
the initial `WAVES$GRAMMSCHMIDT` call into `PAW_GRAM_*` rows. This is a
diagnostic-only change; it does not change the numerical path or runtime
defaults.

Spark C86C validation:

```
runs/gram-profile-final1024-20260531-141705
runs/gram-profile-final2048-20260531-141807
```

| Empty bands | Wall time | `PAW_ETOT_SETUP_GRAM` | `PAW_GRAM_SOLVE` | `PAW_GRAM_PROJECTIONS` | `PAW_GRAM_TRANSFORM` | Final energy |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1024 | 40.38 s | 28.6252 s | 26.5881 s | 0.9343 s | 0.4761 s | 302.280854 Ha |
| 2048 | 279.67 s | 241.9808 s | 237.5816 s | 1.2245 s | 1.5176 s | 302.280854 Ha |

At 2048 bands, the now-fast one-center path is no longer the main problem:
`PAW_1COVERLAP_TOTAL=0.4809 s`, while the initial Gram-Schmidt solve consumes
almost all of `WAVES$ETOT`. The regular orthogonalization phase still has a
smaller solve component (`PAW_ORTHO_SOLVE=8.7036 s` of
`PAW_ORTHO_TOTAL=14.8484 s`), so both call paths point at the same underlying
routine family.

Practical conclusion: the next implementation PR should target
`WAVES_ORTHO_Y_C` or replace the initial Gram-Schmidt orthogonalization with a
more accelerator-friendly dense linear algebra path. Moving more FFT calls to
cuFFT or adding more projector packing will not move the 2048-band wall time
until this solve is addressed.

## Initial Gram-Schmidt Cholesky Solve

The follow-up implementation replaces only the initial
`WAVES$GRAMMSCHMIDT` special case where `PHIPHI=CHIPHI=CHICHI=S`. For a
positive-definite overlap matrix, it computes `S=U^H U` with LAPACK `ZPOTRF`,
uses `ZTRTRI` to form `T=inv(U)`, and applies `X=T-I`. If the Cholesky path
fails, the code falls back to the legacy `WAVES_ORTHO_Y_C` solver. Residency
profile builds enable the path by default; set `CPPAW_GRAM_CHOLESKY=0` to use
the legacy solver for comparison.

Spark C86C validation:

```
runs/gram-cholesky-smoke1024-20260531-143107
runs/gram-cholesky-smoke2048-20260531-143226
runs/gram-cholesky-nstep3-1024-20260531-143332
runs/gram-cholesky-default1024-20260531-143825
runs/gram-cholesky-parallel-smoke512-20260531-144006
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | `PAW_ETOT_SETUP_GRAM` | `PAW_GRAM_SOLVE` | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Legacy default before Cholesky | 1024 | 1 | 1 | 42.10 s | 30.2490 s | 28.2038 s | 302.280854 Ha |
| Cholesky opt-in | 1024 | 1 | 1 | 13.57 s | 2.2011 s | 0.1479 s | 302.280854 Ha |
| New default Cholesky | 1024 | 1 | 1 | 14.05 s | - | - | 302.280854 Ha |
| Explicit legacy override | 1024 | 1 | 1 | 41.74 s | - | - | 302.280854 Ha |
| Legacy default before Cholesky | 2048 | 1 | 1 | 279.67 s | 241.9808 s | 237.5816 s | 302.280854 Ha |
| Cholesky opt-in | 2048 | 1 | 1 | 41.34 s | 5.0733 s | 0.8782 s | 302.280854 Ha |
| Legacy default before Cholesky | 1024 | 3 | 1 | 60.59 s | - | - | 208.886424 Ha |
| Cholesky opt-in | 1024 | 3 | 1 | 32.85 s | - | - | 208.886424 Ha |
| Parallel Cholesky smoke | 512 | 1 | 4 | 9.31 s | - | - | 302.280854 Ha |

The 2048-band result is the clearest design signal: the previous 237.6 s
initial Gram-Schmidt solve shrinks below 0.9 s, and the total run drops from
279.67 s to 41.34 s with the same final energy. This path currently uses CPU
LAPACK through the active NVHPC/NVPL linkage, not a GPU kernel. That is already
enough to remove the dominant bottleneck; a later cuSOLVER `potrf/trtri` variant
is only worth pursuing if larger runs show that this remaining sub-second block
grows again.

## Real Ortho-X Solver Profiling

The next diagnostic patch keeps the numerical path unchanged and splits the
real `WAVES_ORTHO_X` solver inside `PAW_ORTHO_SOLVE` into diagonalization,
residual construction, update, and iteration-count rows. This was added after a
naive exact Cholesky-root experiment was rejected: on the 512-band smoke it
completed but changed the final energy from 302.280854 Ha to 322.671562 Ha.
That means any future exact solve must preserve the physically relevant root
near the current `LAMBDA`, not merely satisfy the quadratic orthogonality
equation.

Spark C86C validation:

```
runs/orthox-profile-smoke512-20260531-145356
runs/orthox-profile-smoke2048-20260531-145404
runs/orthox-profile-parallel-smoke512-20260531-145527
```

| Case | Empty bands | Ranks | Wall time | `PAW_ORTHO_SOLVE` | `PAW_ORTHO_X_DIAG` | `PAW_ORTHO_X_RESIDUAL` | `PAW_ORTHO_X_UPDATE` | Iterations | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial smoke | 512 | 1 | 7.17 s | - | - | - | - | - | 302.280854 Ha |
| Serial large-band smoke | 2048 | 1 | 40.27 s | 8.5702 s | 0.2726 s | 2.8195 s | 5.2874 s | 24 | 302.280854 Ha |
| Parallel smoke | 512 | 4 | 9.31 s | - | - | - | - | - | 302.280854 Ha |

The 2048-band solve now has a precise split: the update half dominates, not the
initial diagonalization. The next implementation attempt should therefore aim
to reduce the 24 Newton-style update iterations or keep the update operands
resident across the transform/back-transform sequence, while preserving the
root selected by the current iterative algorithm.

The follow-up subphase split keeps the same numerical path and divides the
residual and update rows into their matrix products and host-side loops:

```
runs/orthox-update-profile2048-20260531-150019
runs/orthox-update-profile-parallel-smoke512-20260531-150303
```

| Profile row | 2048-band total | Iterations | Interpretation |
| --- | ---: | ---: | --- |
| `PAW_ORTHO_X_DIAG` | 0.2711 s | - | Diagonalization remains small after the Cholesky setup fix. |
| `PAW_ORTHO_X_RESIDUAL` | 2.8238 s | 24 | Residual envelope. |
| `PAW_ORTHO_X_RESIDUAL_MATMUL` | 2.5759 s | 24 | Residual cost is mostly the two BLAS-style products. |
| `PAW_ORTHO_X_RESIDUAL_CHECK` | 0.2479 s | 24 | Symmetrization and convergence check are minor. |
| `PAW_ORTHO_X_UPDATE` | 5.2925 s | 24 | Update envelope. |
| `PAW_ORTHO_X_UPDATE_TRANSFORM` | 2.3523 s | 24 | First half of the update matrix transform. |
| `PAW_ORTHO_X_UPDATE_SCALE` | 0.1937 s | 24 | Eigenvalue-denominator scaling is not the main cost. |
| `PAW_ORTHO_X_UPDATE_BACKTRANSFORM` | 2.3527 s | 24 | Back-transform costs the same as the forward transform. |
| `PAW_ORTHO_X_UPDATE_APPLY` | 0.3671 s | 24 | Host-side `LAMBDA` update is secondary. |
| `PAW_ORTHO_X_UPDATE_SYM` | 0.0266 s | 24 | Occupation symmetrization is negligible. |

The 2048-band run completed in 41.62 s with final energy 302.280854 Ha. The
parallel 512-band smoke completed with 4 MPI ranks in 9.53 s and the same final
energy. This makes the next useful Ortho-X target narrower: optimize or reduce
the repeated transform/back-transform BLAS pairs, not the small scalar loops.

## Ortho-X Workspace Residency

The next implementation prototype adds the opt-in `gpu_resident_orthox` case,
also enabled directly by `CPPAW_GPU_ORTHO_X_RESIDENCY=1`. It keeps
`LAMBDA`, `GAMN`, `HAUX`, `PSIPSI`, `CHIPSI`, `OCC`, and `EIG` in one OpenACC
data region for the real `WAVES_ORTHO_X` iteration loop. The large transform
pairs call the existing present-device cuBLAS helpers directly, while the
residual check, denominator scaling, `LAMBDA` update, and occupation
symmetrization loops run as OpenACC kernels. The numerical update formula is
unchanged.

Spark C86C validation:

```
runs/orthox-resident-smoke512-20260531-151403
runs/orthox-resident-smoke2048-20260531-151428
runs/orthox-resident-parallel-smoke512-20260531-151713
runs/orthox-resident-nstep3-1024-20260531-151840
runs/orthox-resident-ab1024-nstep3-repeat3-20260531-152238
runs/orthox-resident-ab2048-nstep3-20260531-152609
runs/orthox-cublas-accounting-smoke512-final-20260531-153740
runs/orthox-cublas-accounting-parallel-smoke512-20260531-153831
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `PAW_ORTHO_X_RESIDUAL` | `PAW_ORTHO_X_UPDATE` | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.15 s | 3.6863 GB | - | - | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.70 s | 3.0571 GB | - | - | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 41.63 s | 22.1456 GB | 2.8577 s | 5.3654 s | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 39.86 s | 8.2058 GB | 2.4293 s | 4.6135 s | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 4 | 9.38 s | 3.5241 GB | - | - | 302.280854 Ha |
| `gpu_resident_orthox` | 1024, `NSTEPS=3` | 1 | 33.42 s | 9.8012 GB | - | - | 208.886424 Ha |
| `gpu_resident`, repeat avg | 1024, `NSTEPS=3` | 1 | 33.99 s | 18.8042 GB | - | - | 208.886424 Ha |
| `gpu_resident_orthox`, repeat avg | 1024, `NSTEPS=3` | 1 | 32.79 s | 9.8012 GB | - | - | 208.886424 Ha |
| `gpu_resident` | 2048, `NSTEPS=3` | 1 | 105.64 s | 56.2853 GB | - | - | 208.886424 Ha |
| `gpu_resident_orthox` | 2048, `NSTEPS=3` | 1 | 106.04 s | 17.4964 GB | - | - | 208.886424 Ha |

The 2048-band case is the useful signal: wall time improves by about 4.3%, the
copy estimate drops by about 13.9 GB, the residual check shrinks from 0.2489 s
to 0.0207 s, and the update scalar/apply/sym loops shrink from about 0.589 s to
0.0188 s. The two transform/back-transform matrix-product pairs remain the
dominant Ortho-X work. The accounting follow-up records the direct
present-device cuBLAS calls as `CUBLAS_DGEMM_ORTHOX_RESIDUAL`,
`CUBLAS_DGEMM_ORTHOX_TRANSFORM`, and `CUBLAS_DGEMM_ORTHOX_BACKTRANS`, so the
benchmark `blas_s` column now includes this prototype. The accounting smoke
checks kept the 512-band energies unchanged at 302.280854 Ha for both serial
and 4-rank runs. Keep the path opt-in until longer runs check multi-step energy
stability and larger systems.

The longer follow-up keeps that conclusion nuanced. At 1024 empty bands and
`NSTEPS=3`, three repeats show a consistent wall-time gain for
`gpu_resident_orthox` (32.79 s average versus 33.99 s). At 2048 empty bands and
`NSTEPS=3`, the same path is wall-time neutral on Spark (106.04 s versus
105.64 s) while still cutting the copy estimate from 56.3 GB to 17.5 GB. The
standard NVHPC comparison now includes `gpu_resident_orthox` so future runs keep
tracking this tradeoff, but the recommended default remains plain
`gpu_resident`.

### Nsight Case Harness

The Nsight Systems harness now accepts a `CASE` keyword, selects the matching
serial or parallel executable, applies the GPU/library runtime switches, and
writes the resolved settings to `nsys_case.env` in the run directory. Spark
C86C smoke traces:

```
runs/nsys-case-gpu_resident-smoke512-20260531-154536
runs/nsys-case-gpu_resident_orthox-smoke512-20260531-154549
runs/nsys-case-gpu_resident_orthox-final-smoke512-20260531-154929
```

Both 512-band traces produced `nsys.nsys-rep`, `nsys_profile.csv`, and
`summary.txt`. The `gpu_resident_orthox` trace recorded
`CPPAW_GPU_ORTHO_X_RESIDENCY=1` in `nsys_case.env`, kept the energy at
302.280854 Ha, and included `CUBLAS_DGEMM_ORTHOX_RESIDUAL`,
`CUBLAS_DGEMM_ORTHOX_TRANSFORM`, and `CUBLAS_DGEMM_ORTHOX_BACKTRANS` in the
profile CSV.

The 2048-band Nsight comparison made the synchronization tradeoff explicit:

```
runs/nsys-case-gpu_resident-smoke2048-20260531-155236
runs/nsys-case-gpu_resident_orthox-smoke2048-20260531-155323
runs/nsys-case-gpu_resident_orthox_nosync-smoke2048-20260531-155936
runs/orthox-nosync-case2048-20260531-155722
```

| Case | Wall time | Instrumented rank-s | Gap | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_orthox` | 39.54 s | 24.7906 s | 14.7494 s | 8.2058 GB | 302.280854 Ha |
| `gpu_resident_orthox_nosync` | 39.73 s | 17.2575 s | 22.4725 s | 8.2058 GB | 302.280854 Ha |
| `gpu_resident_nosync` | 41.14 s | 24.4806 s | 16.6594 s | 22.1456 GB | 302.280854 Ha |

`gpu_resident_orthox_nosync` removes the explicit cuBLAS-side
`cudaDeviceSynchronize` time from the instrumented rows, but the Nsight trace
then shows the waiting time reappearing in `cuStreamSynchronize` and
device-to-host runtime calls. The case is therefore useful for profiling
attribution, but it is not a new recommended default on Spark.

The Nsight harness also writes a compact SQL-derived CUDA activity summary when
Nsight leaves a `.sqlite` export. Smoke run
`runs/nsys-sql-summary-smoke512-20260531-160604` produced
`nsys_sql_summary.txt` with top kernel, runtime, memcpy, and synchronization
tables; for the `gpu_resident_orthox_nosync` smoke it highlighted
`cuStreamSynchronize` as the leading runtime row and kept the energy unchanged
at 302.280854 Ha.

### Present-Aware cuBLAS Copy Accounting

The generic resident cuBLAS wrappers now split estimated transfers into
per-array `ACC_PRESENT_CUBLAS_*` and `ACC_COPY_CUBLAS_*` rows instead of using
the older aggregate `ACC_COPY_CUBLAS_{ZSPROD,DSPROD,ZGEMM_NN}_RES` estimates.
For Hermitian/symmetric scalarproduct calls with `TID=.TRUE.`, the OpenACC data
region also omits the unused second wavefunction array. This is primarily a
profiling correctness change: it keeps already-resident arrays visible without
charging them as host/device copies, and it leaves the numerical path
unchanged.

Spark C86C validation:

```
runs/present-aware-cublas-copy-smoke512-20260531-161837
runs/present-aware-cublas-copy-parallel512-20260531-161936
runs/present-aware-cublas-copy2048-20260531-162011
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Generic resident cuBLAS copy rows | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.25 s | 4.1716 GB | 1.4993 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.09 s | 3.4606 GB | 0.7103 GB | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.40 s | 6.7718 GB | 0.5943-0.5949 GB per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 41.71 s | 25.2977 GB | 10.0257 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 40.05 s | 9.5775 GB | 4.6847 GB | 302.280854 Ha |

All checked profile CSV files contained the new granular rows and no hits for
the old aggregate `ACC_COPY_CUBLAS_ZSPROD_RES`,
`ACC_COPY_CUBLAS_DSPROD_RES`, or `ACC_COPY_CUBLAS_ZGEMM_NN_RES` rows. In the
2048-band orthox case, the remaining large generic cuBLAS copy estimates are
now explicit: `ACC_COPY_CUBLAS_ZGEMM_NN_C_IO` accounts for 2.7434 GB, while the
real `DSPROD` input/output rows account for the next major block. That makes the
next residency target clearer than the old aggregate estimate did.

### Residency Copy Accounting Cleanup

A follow-up removed the older aggregate region estimates
`ACC_COPY_CUBLAS_PROJ_RES`, `ACC_COPY_CUBLAS_OVERLAP_RES_REGION`, and
`ACC_COPY_CUBLAS_ADDOPSI_RES_REGION`. Their component arrays are already covered
by the present-aware `ACC_PRESENT_*` / `ACC_COPY_*` rows, so keeping the
aggregate rows double-counted the same residency regions in the `copy_gb`
summary. The kernel path is unchanged; this only makes the copy accounting match
the array-level profile rows.

Spark C86C validation:

```
runs/residency-copy-accounting-cleanup512-20260531-162946
runs/residency-copy-accounting-cleanup2048-20260531-163012
runs/residency-copy-accounting-cleanup-parallel512-20260531-163144
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.49 s | 2.9626 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.98 s | 2.2516 GB | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.42 s | 5.4861 GB | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 42.15 s | 23.7046 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 38.94 s | 7.9844 GB | 302.280854 Ha |

All checked CSV files had no hits for the removed aggregate rows. After cleanup,
the largest remaining explicit copy target in the 2048-band orthox profile is
`ACC_COPY_CUBLAS_ZGEMM_NN_C_IO` at 2.7434 GB, followed by real wavefunction
input/output rows such as `ACC_COPY_PROJ_PSI_IN`, `ACC_COPY_GRAM_PSI_IN`, and
`ACC_COPY_ADDPRO_PSI_IO` at 0.4572 GB each.

### Gram Transform Residency

The Gram-Schmidt resident region now keeps `PSI` on the device through the final
wavefunction transform, changing the outer Gram row from `ACC_COPY_GRAM_PSI_IN`
to `ACC_COPY_GRAM_PSI_IO`. This lets the `LIB$ADDPRODUCTC8` transform calls see
their output matrix as present instead of opening an additional cuBLAS copy
region for the same wavefunction data.

Spark C86C validation:

```
runs/gram-transform-residency512-20260531-163751
runs/gram-transform-residency2048-20260531-163818
runs/gram-transform-residency-parallel512-20260531-163945
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `ZGEMM_NN_C_IO` | `ZGEMM_NN_C` present | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.14 s | 2.5592 GB | 2 calls, 0.2690 GB | 4 calls | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.24 s | 1.8481 GB | 2 calls, 0.2690 GB | 4 calls | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.31 s | 5.0827 GB | 2 calls, 0.0672 GB per rank | 4 calls per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 40.26 s | 22.3328 GB | 2 calls, 0.9145 GB | 4 calls | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 39.13 s | 6.6127 GB | 2 calls, 0.9145 GB | 4 calls | 302.280854 Ha |

Compared with the previous cleanup run, the 2048-band orthox copy estimate drops
from 7.9844 GB to 6.6127 GB. The remaining `ZGEMM_NN_C_IO` calls are therefore
not Gram-transform output copies; they are the next addproduct-residency target.

### ADDOPSI Output Residency

The orthogonalization `WAVES_ADDOPSI` step now keeps the updated wavefunction
array resident while `LIB$ADDPRODUCTC8` accumulates into it. This removes the
last `ACC_COPY_CUBLAS_ZGEMM_NN_C_IO` rows from the Si64 inversion-symmetric
path. To avoid stale device data, the inversion path keeps only the output
`PSIM`/`PSIBAR` array in the OpenACC data region; `OPSI` and the temporary
lambda blocks are still copied by the per-call cuBLAS wrapper because `OPSI` is
inverted on the host between the two addproduct calls.

Spark C86C validation:

```
runs/addopsi-output-residency512-20260531-164747
runs/addopsi-output-residency2048-20260531-164801
runs/addopsi-output-residency-parallel512-20260531-164922
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `ADDOPSI_PSIM_IO` | `ZGEMM_NN_C_IO` | `ZGEMM_NN_C` present | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 6.66 s | 2.4247 GB | 1 call, 0.1345 GB | 0 calls | 6 calls | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.17 s | 1.7136 GB | 1 call, 0.1345 GB | 0 calls | 6 calls | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.32 s | 4.9482 GB | 1 call, 0.0336 GB per rank | 0 calls | 6 calls per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 41.63 s | 21.8756 GB | 1 call, 0.4572 GB | 0 calls | 6 calls | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 39.04 s | 6.1555 GB | 1 call, 0.4572 GB | 0 calls | 6 calls | 302.280854 Ha |

Compared with the Gram-transform residency run, the 2048-band orthox copy
estimate drops from 6.6127 GB to 6.1555 GB. The total drop is smaller than the
old `ZGEMM_NN_C_IO` estimate because the remaining required output transfer is
now charged explicitly to `ACC_COPY_ADDOPSI_PSIM_IO`.

### Gram Transform Device Input

The final Gram-Schmidt transform now creates the `PSIINV` scratch wavefunction
on the device and fills it from the already resident `PSI` array before calling
`PLANEWAVE$ADDPRODUCT`. The host `PSIINV=PSI` assignment is kept for CPU
fallback correctness, but the residency-profile cuBLAS path no longer needs a
host/device copy for this scratch input. The transform matrices `X`, `X1`, and
`X2` are still host-computed and copied explicitly.

Spark C86C validation:

```
runs/gram-transform-device-input512-20260531-170001
runs/gram-transform-device-input2048-20260531-170014
runs/gram-transform-device-input-parallel512-20260531-170136
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `PSIINV` row | Remaining `ZGEMM_NN_A_IN` | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 6.39 s | 2.2902 GB | present, 2 calls | 2 calls, 0.1345 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.05 s | 1.5792 GB | present, 2 calls | 2 calls, 0.1345 GB | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.45 s | 4.8137 GB | present, 2 calls per rank | 2 calls, 0.0336 GB per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 41.23 s | 21.4184 GB | present, 2 calls | 2 calls, 0.4572 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 40.20 s | 5.6982 GB | present, 2 calls | 2 calls, 0.4572 GB | 302.280854 Ha |

Compared with the ADDOPSI output residency run, the 2048-band orthox copy
estimate drops from 6.1555 GB to 5.6982 GB. The remaining
`ACC_COPY_CUBLAS_ZGEMM_NN_A_IN` rows now correspond to the `WAVES_ADDOPSI`
`OPSI` inputs, not the Gram transform scratch copy.

### ADDOPSI Device Inversion

The inversion-symmetric `WAVES_ADDOPSI` path now uses the same inversion-batch
idea as the plane-wave addproduct helper. The residency-profile GPU path copies
`OPSI` once, builds an inverted `OPSIINV` scratch array on the device using
`MINUSG`, and runs both addproducts with present inputs. The original host
inversion branch remains the fallback when residency or inversion batching is
disabled.

Spark C86C validation:

```
runs/addopsi-device-inversion512-20260531-170700
runs/addopsi-device-inversion2048-20260531-170715
runs/addopsi-device-inversion-parallel512-20260531-170835
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `OPSI_TINV_IN` | `OPSIINV_TINV` | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.27 s | 2.2230 GB | 1 call, 0.0672 GB | present, 1 call | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.41 s | 1.5119 GB | 1 call, 0.0672 GB | present, 1 call | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.32 s | 4.7465 GB | 1 call, 0.0168 GB per rank | present, 1 call per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 40.98 s | 21.1897 GB | 1 call, 0.2286 GB | present, 1 call | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 39.37 s | 5.4696 GB | 1 call, 0.2286 GB | present, 1 call | 302.280854 Ha |

Compared with the Gram-transform device-input run, the 2048-band orthox copy
estimate drops from 5.6982 GB to 5.4696 GB. The remaining large wavefunction
copy rows are now outside this local addproduct path: the outer Gram
`PSI` input/output region, projection wavefunction inputs, and force/setup
wavefunction copies.

### Orthogonalization Wave Region

The main orthogonalization OpenACC data region now keeps `PSIM` and `OPSI`
resident from the projection/overlap phase through the final `WAVES_ADDOPSI`
update. `PSIM` is copied back once at the end of that broader region, and
`OPSI` remains input-only. This removes the separate ADDOPSI output/input
copies while keeping the CPU Lagrange solve unchanged.

Spark C86C validation:

```
runs/orthogonalize-wave-region512-20260531-171440
runs/orthogonalize-wave-region2048-20260531-171455
runs/orthogonalize-wave-region-parallel512-20260531-171612
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | ADDOPSI wave rows | Final energy |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `gpu_resident` | 512 | 1 | 7.66 s | 2.0885 GB | `PSIM` present, `OPSI` present | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 6.77 s | 1.3774 GB | `PSIM` present, `OPSI` present | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.25 s | 4.6120 GB | `PSIM` present, `OPSI` present per rank | 302.280854 Ha |
| `gpu_resident` | 2048 | 1 | 39.02 s | 20.7325 GB | `PSIM` present, `OPSI` present | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 38.04 s | 5.0124 GB | `PSIM` present, `OPSI` present | 302.280854 Ha |

Compared with the ADDOPSI device-inversion run, the 2048-band orthox copy
estimate drops from 5.4696 GB to 5.0124 GB. `ACC_COPY_ADDOPSI_PSIM_IO` and
`ACC_COPY_ADDOPSI_OPSI_TINV_IN` disappear from the checked profiles; the broader
region instead records one `ACC_COPY_ORTHO_PSIM_IO` and one
`ACC_COPY_ORTHO_OPSI_IN`.

### ADDPRO Context Copy Accounting

`WAVES_ADDPRO` now tags its residency profile rows by caller context. The old
aggregate rows `ACC_COPY_ADDPRO_PSI_IO`, `ACC_COPY_ADDPRO_PROPSI_IN`, and
`ACC_PRESENT_ADDPRO_PRO_CACHE` are split into `HPSI` and `OPSI` rows, for
example `ACC_COPY_ADDPRO_HPSI_PSI_IO` and `ACC_COPY_ADDPRO_OPSI_PSI_IO`.
The `PSI` row also uses the input/output present-check helper, matching the
`PRESENT_OR_COPY(PSI)` OpenACC data region used by the projector addition.

Spark C86C validation:

```
runs/addpro-profile-contexts512-20260531-172542
runs/addpro-profile-contexts-parallel512-20260531-172641
runs/addpro-profile-contexts2048-20260531-172651
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | ADDPRO `PSI_IO` rows | Final energy |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `gpu_resident` | 512 | 1 | 6.66 s | 2.2230 GB | `HPSI` 0.1345 GB, `OPSI` 0.1345 GB | 302.280854 Ha |
| `gpu_resident_orthox` | 512 | 1 | 7.34 s | 1.5119 GB | `HPSI` 0.1345 GB, `OPSI` 0.1345 GB | 302.280854 Ha |
| `gpu_resident` | 512 | 4 | 9.52 s | 4.7465 GB | context rows per rank | 302.280854 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 38.88 s | 5.4696 GB | `HPSI` 0.4572 GB, `OPSI` 0.4572 GB | 302.280854 Ha |

This is an accounting and localization change, not a new optimization. The
2048-band orthox copy estimate rises from 5.0124 GB to 5.4696 GB because
`ADDPRO_PSI` is now counted as copy-in plus copy-out. The split shows the next
real residency target clearly: one updated wavefunction copy comes from the
Hamiltonian application (`HPSI`) and one from overlap-wave construction (`OPSI`).

## Energy-Guard Default Refresh

After the Si64 energy guard was added, Spark C86C was rebuilt from
`cp-paw-nvhpc` and rerun with the current residency implementation:

```
runs/fresh-energy-guard-sweep-20260531-174614
```

Focused `TEST=si64`, `NSTEPS=1` smokes:

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.12 s | 2.2230 GB | 0.000000401 Ha |
| `gpu_resident_orthox` | 512 | 1 | 6.90 s | 1.5119 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.31 s | 4.7465 GB | 0.000000401 Ha |
| `gpu_resident_orthox` | 512 | 4 | 9.19 s | 1.9022 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 41.92 s | 21.1897 GB | 0.000000407 Ha |
| `gpu_resident_orthox` | 2048 | 1 | 40.89 s | 5.4696 GB | 0.000000407 Ha |

The refreshed `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1` standard
comparison gives:

| Suite | Case | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | --- | ---: | ---: | ---: | ---: |
| GPU | `gpu_resident_orthox` | 1 | 13.39 s | 2.7137 GB | 0.000000407 Ha |
| GPU | `gpu_resident_addpro_host` | 1 | 13.71 s | 6.2208 GB | 0.000000407 Ha |
| GPU | `gpu_resident` | 1 | 13.74 s | 6.3553 GB | 0.000000407 Ha |
| GPU | `gpu_resident_pro_host` | 1 | 14.16 s | 7.2698 GB | 0.000000407 Ha |
| GPU | `gpu_resident_no_cusolver` | 1 | 17.02 s | 6.2597 GB | 0.000000407 Ha |
| GPU | `gpu_off` | 1 | 74.59 s | 0.0000 GB | 0.000000398 Ha |
| CPU | `cpu` | 1 | 70.83 s | 0.0000 GB | 0.000000398 Ha |
| CPU | `nvhpc_cpu` | 1 | 71.85 s | 0.0000 GB | 0.000000398 Ha |
| CPU ref | `cpu` | 8 | 168.36 s | 0.0000 GB | 0.000000398 Ha |
| CPU ref | `nvhpc_cpu` | 8 | 167.20 s | 0.0000 GB | 0.000000398 Ha |

These runs keep the energy check green and make the default change less
speculative: `gpu_resident_orthox` is both copy-lighter and slightly faster at
512, 2048 and the 1024-band standard point. The residency-profile build now
enables `CPPAW_GPU_ORTHO_X_RESIDENCY=1` by default; use
`gpu_resident_orthox_off` or `CPPAW_GPU_ORTHO_X_RESIDENCY=0` to compare against
the previous path.

Patch validation for the default flip:

```
runs/orthox-default-smoke-20260531-180220
```

| Case | Wall time | Total copy estimate | Energy delta | Meaning |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident` | 7.22 s | 1.5119 GB | 0.000000401 Ha | Default now takes the ORTHO_X resident path. |
| `gpu_resident_orthox_off` | 7.24 s | 2.2230 GB | 0.000000401 Ha | Explicit comparison with the previous path. |

## Gram Context Copy Accounting

`WAVES_GRAMSCHMIDT` now tags its outer resident wavefunction copy by caller
context. The old aggregate `ACC_COPY_GRAM_PSI_IO` row is split into
`ACC_COPY_GRAM_PSI0_PSI_IO` and `ACC_COPY_GRAM_PSIM_PSI_IO`, while the scratch
and transform rows remain unchanged.

Spark C86C validation:

```
runs/gram-profile-contexts-20260531-180922
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Gram `PSI_IO` rows | Energy delta |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `gpu_resident` | 512 | 1 | 6.94 s | 1.5119 GB | `PSI0` 0.1345 GB, `PSIM` 0.1345 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.15 s | 1.9022 GB | context rows per rank | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 41.02 s | 5.4696 GB | `PSI0` 0.4572 GB, `PSIM` 0.4572 GB | 0.000000407 Ha |

This is instrumentation, not a speedup. It shows that the remaining initial
Gram wavefunction transfer is symmetric between the `PSI0` and `PSIM`
orthogonalization calls. A later optimization should therefore either cover both
initial Gram calls with a broader, verified wavefunction lifetime or leave this
path alone and focus first on the non-Gram projection/Addpro edges.

## Projection Context Copy Accounting

`WAVES_PROJECTIONS` now tags the profiled `PSI` input copy by caller context.
The old aggregate `ACC_COPY_PROJ_PSI_IN` row is split into rows such as
`ACC_COPY_PROJ_SETUP0_PSI_IN`, `ACC_COPY_PROJ_WRITEPDOS_PSI_IN`, and present
checks for resident callers such as `ACC_PRESENT_PROJ_GRAM_PSI0_PSI`.
`PROPSI_OUT` stays aggregate because the profiler label length is intentionally
kept short.

Spark C86C validation:

```
runs/projection-profile-contexts-20260531-181555
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Projection `PSI` rows | Energy delta |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `gpu_resident` | 512 | 1 | 7.04 s | 1.5119 GB | `SETUP0` 0.0672 GB, `WRITEPDOS` 0.0672 GB; Gram/Ortho present | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.36 s | 1.9022 GB | `SETUP0` and `WRITEPDOS` about 0.0168 GB per rank; Gram/Ortho present | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 39.94 s | 5.4696 GB | `SETUP0` 0.2286 GB, `WRITEPDOS` 0.2286 GB; Gram/Ortho present | 0.000000407 Ha |

This is instrumentation, not a speedup. The split moves an important design
decision out of the dark: the remaining projection wavefunction copies in the
current Si64 smoke are setup and PDOS/reporting edges, while the Gram and
orthogonalization projection calls already see resident `PSI`. The next
optimization should therefore prioritize the measured `ADDPRO_HPSI` and
`ADDPRO_OPSI` edges before widening projection residency.

## Overlap Context Copy Accounting

`WAVES_OVERLAP` now tags the generic cuBLAS `ZSPROD` copy/present rows by
caller context. The old aggregate rows such as
`ACC_COPY_CUBLAS_ZSPROD_PSI1_IN`, `ACC_COPY_CUBLAS_ZSPROD_PSI2_IN`,
`ACC_COPY_CUBLAS_ZSPROD_OUT`, and `ACC_COPY_CUBLAS_ZSPROD_OVL_RES` are split
into short rows such as `ACC_COPY_ZSP_HAMILTON_P1_IN`,
`ACC_PRESENT_ZSP_GRAM_PSI0_P1`, and `ACC_COPY_ZSP_ORTH_PSIM_OVL`.

Spark C86C validation:

```
runs/overlap-profile-contexts-20260531-162835
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.03 s | 1.5119 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.35 s | 1.9022 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 39.04 s | 5.4696 GB | 0.000000407 Ha |

The 2048-band context split shows that the large remaining non-resident overlap
copies are Hamiltonian/reporting edges:

| Row | Calls | Copy estimate | Meaning |
| --- | ---: | ---: | --- |
| `ACC_COPY_ZSP_HAMILTON_P1_IN` | 2 | 0.4572 GB | First Hamiltonian overlap input is copied. |
| `ACC_COPY_ZSP_HAMILTON_P2_IN` | 1 | 0.2286 GB | Second Hamiltonian overlap input is copied once; the other call sees it present. |
| `ACC_COPY_ZSP_HAMILTON_OUT` | 2 | 0.0379 GB | Hamiltonian overlap output copy. |
| `ACC_PRESENT_ZSP_GRAM_PSI0_P1/P2` | 1 each | 0 GB | Gram overlap inputs are resident. |
| `ACC_PRESENT_ZSP_GRAM_PSIM_P1/P2` | 1 each | 0 GB | Gram overlap inputs are resident. |
| `ACC_PRESENT_ZSP_ORTH_*_P1/P2` | 1 each | 0 GB | Orthogonalization overlap inputs are resident. |

This is instrumentation, not a speedup. It shows that widening overlap
residency inside Gram/Ortho would not attack the dominant remaining `ZSPROD`
copies in this smoke. The next optimization target remains the measured
wavefunction updates in `ACC_COPY_ADDPRO_HPSI_PSI_IO` and
`ACC_COPY_ADDPRO_OPSI_PSI_IO`; Hamiltonian/reporting residency can be considered
later if the band-output path matters for production runs.

## OPSI Build-Residency Diagnostic

`CPPAW_GPU_OPSI_RESIDENCY=1` now enables an opt-in path that can keep
orthogonalization `OPSI` resident from its `WAVES_OPSI` build through the later
projection, overlap, and `WAVES_ADDOPSI` phase when the path is non-stress,
non-superwave, and all atom blocks pass the addproduct threshold. The Si64 smoke
uses inversion-symmetric superwaves (`NBH != NB`), so the correctness guard
deliberately leaves the committed diagnostic inactive for this case.

Spark C86C validation:

```
runs/opsi-build-residency-20260531-164302
runs/opsi-build-residency-20260531-165523  # exact final-commit 512/1 check
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.82 s | 1.5119 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 1 | 7.64 s | 1.5119 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.30 s | 1.9022 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 4 | 9.19 s | 1.9022 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 40.92 s | 5.4696 GB | 0.000000407 Ha |
| `gpu_resident_opsi` | 2048 | 1 | 40.94 s | 5.4696 GB | 0.000000407 Ha |

The exact final commit was also rebuilt and checked at 512/1:
`gpu_resident` 6.99 s versus `gpu_resident_opsi` 7.25 s, both with
1.5119 GB copy estimate and a valid 0.000000401 Ha energy delta.

A discarded superwave-widening experiment did reduce the copy estimate
(`gpu_resident_opsi` at 512/1: 1.3774 GB; at 2048/1: 5.0124 GB), but it changed
the Si64 energy to 296.752801 Ha, about 5.528 Ha away from the reference. That
made the resident superwave overlap path the next correctness target, followed
by a narrower retest of where superwave OPSI can safely enter the resident
region.

## Superwave Overlap Residency

The follow-up patch moves the superwave `<PSI_-|PSI_+>` inversion contribution
onto an explicit resident cuBLAS path, recorded as `CUBLAS_ZGEMM_OVL_RES_INV`.
It also makes the gamma correction in the resident `<PSI_+|PSI_+>` branch read
the gamma slice from device data. This removes the last host-side overlap
assumption that blocked a safe superwave OPSI staging experiment.

The fully resident superwave OPSI build/scale path was still energy-invalid:
building OPSI with resident `WAVES_ADDPRO`, or scaling the freshly built OPSI on
the device, reproduced the bad 296.752801 Ha Si64 energy. The committed
superwave OPSI diagnostic is therefore conservative: build OPSI on the host,
mass-scale it on the host, and only then enter the device-resident
projection/overlap/`WAVES_ADDOPSI` region. This is correctness-preserving but
does not reduce the Si64 copy estimate versus `gpu_resident`; it is useful as a
guarded staging point for future superwave `WAVES_ADDPRO` and mass-scaling work.

Spark C86C validation:

```
runs/superwave-opsi-hostscale-20260531-512-1r
runs/superwave-opsi-hostscale-20260531-512-4r
runs/superwave-opsi-hostscale-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.09 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 1 | 7.06 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.45 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 4 | 9.15 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 38.74 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_opsi` | 2048 | 1 | 38.73 s | 5.3749 GB | 0.000000407 Ha |

## ADDPRO Context-Cache Controls

The next diagnostic split adds independent HPSI and OPSI controls for
`WAVES_ADDPRO` cache reuse:

- `CPPAW_GPU_ADDPRO_CACHE_HPSI=0`
- `CPPAW_GPU_ADDPRO_CACHE_OPSI=0`

This lets the superwave OPSI follow-up isolate whether a future bad energy comes
from the Hamiltonian-side projector update or from the overlap-wave OPSI update,
without disabling the shared resident `PRO` cache for projections.

Spark C86C validation:

```
runs/addpro-context-cache-20260531-512-1r
runs/addpro-context-cache-20260531-512-4r
runs/addpro-context-cache-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.46 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_addpro_hpsi_host` | 512 | 1 | 6.35 s | 1.5441 GB | 0.000000401 Ha |
| `gpu_resident_addpro_opsi_host` | 512 | 1 | 6.37 s | 1.5441 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 1 | 6.41 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_opsi_addpro_host` | 512 | 1 | 7.24 s | 1.5441 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.56 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_addpro_hpsi_host` | 512 | 4 | 9.30 s | 1.9098 GB | 0.000000401 Ha |
| `gpu_resident_addpro_opsi_host` | 512 | 4 | 9.34 s | 1.9098 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 4 | 9.47 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_opsi_addpro_host` | 512 | 4 | 9.37 s | 1.9098 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 40.42 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_addpro_hpsi_host` | 2048 | 1 | 40.46 s | 5.0925 GB | 0.000000407 Ha |
| `gpu_resident_addpro_opsi_host` | 2048 | 1 | 40.49 s | 5.0925 GB | 0.000000407 Ha |
| `gpu_resident_opsi` | 2048 | 1 | 40.03 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_opsi_addpro_host` | 2048 | 1 | 40.27 s | 5.0925 GB | 0.000000407 Ha |

The switches do not change the Si64 energy. Wall time is essentially neutral
within one-run noise. The 2048 copy estimate decreases when the cached ADDPRO
context is disabled because the current profiler counts the cached path's
`PSI_IO` residency edge explicitly, while the fallback path is less granular.
That makes these switches useful for isolation, not a reason to change the
recommended default.

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
  GPU_CASES="gpu_off gpu_all gpu_all_off gpu_resident gpu_resident_orthox_off gpu_resident_no_cusolver cublas cusolver cufft cufftw nvlamath nvblas gpu_no_cufft gpu_no_cublas gpu_no_cusolver gpu_managed gpu_unified" \
  ./run_nvhpc_standard.sh
```
