# NVHPC Spark Si64 Benchmark Summary

This note summarizes the Spark C86C Si64 band benchmarks used to choose the
current NVIDIA HPC SDK development defaults. Unless noted otherwise, the case
is the periodic Si64 profile input with `TEST=si64_bands`,
`EMPTY_BANDS=1024`, `NSTEPS=1`, one GPU rank for GPU cases, and one-rank plus
eight-rank CPU references.

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
| `1cov-split-20260531-*` | One-center overlap copy accounting | `gpu_resident_hpsi` 40.48 s at 2048/1 | - | `gpu_resident_hpsi` 9.78 s at 512/4 | Splits the 1COV copy estimate into packed projector input and overlap-matrix output rows. |
| `wave-io-split-20260531-*` | Wavefunction IO copy accounting | `gpu_resident_hpsi` 39.22 s at 2048/1 | - | `gpu_resident_hpsi` 9.72 s at 512/4 | Splits Gram, ORTHO, ADDPRO, and ADDOPSI wavefunction IO estimates into input and output rows. |
| `denmat-energy-acc-v2-20260531-*` | DENMAT energy OpenACC diagnostic | `gpu_resident_hpsi` 39.00 s, `gpu_resident_hpsi_denmat_energy` 39.40 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy` 9.75 s at 512/4 | Adds an opt-in two-stage OpenACC diagnostic for the time-inversion DENMAT energy/Lambda contraction; DENMAT shrinks, but Lambda copies keep it diagnostic-only. |
| `denmat-lagr-residency-20260531-*` | DENMAT LAGR setup/residency | `gpu_resident_hpsi` 37.19 s, `gpu_resident_hpsi_denmat_energy` 37.40 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy` 9.56 s at 512/4 | Precomputes `LAGR=LAMBDA*OCC` once per k-point/spin and keeps it resident for the DENMAT diagnostic; copy volume drops sharply, wall time remains neutral at 2048/1. |
| `offden-profile-split-20260531-*` | Off-site DENMAT profiling split | `gpu_resident_hpsi_denmat_energy` 36.39 s at 2048/1 | - | `gpu_resident_hpsi` 9.75 s at 512/4 | Splits `PAW_OFFDEN_SUM` into setup/zero/local/combine; off-site time is mostly local contraction, not MPI combine. |
| `offden-blas-diagnostic-20260531-*` / `offden-blas-combined-20260531-*` | Off-site DENMAT BLAS diagnostic | `gpu_resident_hpsi_denmat_energy_offden_blas` 36.91 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_blas` 9.76 s at 512/4 | Rewrites scalar `TINV` off-site local work as packed `ZGEMM`; energy-valid, much faster inside `PAW_OFFDEN_SUM_LOCAL`, still opt-in host-data diagnostic. |
| `offden-cublas-diagnostic-20260531-*` / `offden-cublas-combined-20260531-*` | Naive off-site DENMAT cuBLAS diagnostic | `gpu_resident_hpsi_denmat_energy_offden_blas` remains better at 36.79 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_blas` remains better at 9.75 s at 512/4 | Adds a forced per-neighbor cuBLAS diagnostic; energy-valid, but kernel timings are worse than host BLAS, so the next GPU attempt must batch or keep buffers resident. |
| `offden-cublas-stack-diagnostic-20260531-*` / `offden-cublas-stack-combined-20260531-*` | Stacked off-site DENMAT cuBLAS diagnostic | `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` 36.10 s at 2048/1 | - | `gpu_resident_hpsi_offden_cublas_batch` 9.65 s at 512/4 | Groups neighbors with the same first atom and second-projector size into one wider cuBLAS `ZGEMM`; energy-valid and eliminates the tiny-GEMM launch problem, but host packing still dominates enough to keep it opt-in. |
| `offden-device-pack-20260531-*` / `offden-device-pack-combined-20260531-*` | Off-site DENMAT device-pack diagnostic | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` 36.03 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` 9.87 s at 512/4 | Packs the stacked off-site A/B buffers on the GPU and copies back only WORK; strong for 1 MPI/GPU, diagnostic-only when several ranks share one GPU. |
| `proj-residency-fixed-20260531-*` | `THIS%PROJ` residency diagnostic | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` 36.57 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` 9.48 s at 512/4 | Keeps the combined projection result present for eligible off-site device-pack consumers; energy-valid after invalidating stale present `PROPSI`, useful as an opt-in diagnostic but too narrow for default promotion. |
| `offden-device-accum-20260601-*` | Off-site DENMAT device-accum diagnostic | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` 35.26 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` 9.59 s at 512/4 | Converts stacked complex `WORK` to real `MATPACK` on the GPU and copies that back; energy-valid and reduces copy volume, but kernel overhead keeps it opt-in. |
| `offden-flat-accum-20260601-*` | Off-site DENMAT flat device-accum diagnostic | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` 35.64 s at 2048/1 | - | `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum` 9.62 s at 512/4 | Accumulates batches into one flat real off-site matrix buffer on the GPU and copies it back once per k-point/spin pass; energy-valid, useful at 1 MPI/GPU, neutral/noisy when four ranks share one GPU. |
| `psim-present-copy-20260601-*` | PSIM propagation present-or-copy accounting | `gpu_resident` 6.36 s at 512/1 | - | `gpu_resident` 9.18 s at 512/4 | Changes the PSIM propagation data region to `present_or_copy`; energy-valid and correct for future broader residency, but current Si64 still copies PSIM in/out because no enclosing resident producer is active. |
| `psim-focus-harness-20260601-*` | Focused PSIM propagation harness | `gpu_psim_propagate` 6.20 s at 512/1 | - | `gpu_psim_propagate` 9.19 s at 512/4 | Adds a reusable PSIM/HPSI propagation sweep; all cases are energy-valid, but the single-run 512-band timings are noisy and the PSIM path still increases copy volume, so this is a regression harness rather than a default promotion. |
| `psim-phase-residency-20260601-*` | Cross-phase PSIM residency diagnostic | `gpu_resident` 37.61 s, `gpu_resident_psim_phase` 37.75 s at 2048/1 | - | `gpu_resident_psim_phase` 9.11 s at 512/4 | Leaves propagated `PSIM` resident into orthogonalization and copies it back at the orthogonalization boundary; energy-valid and reduces copy volume, but wall time is neutral/noisy, so keep it opt-in. |
| `psim-lifecycle-20260601-rerun-*` | Two-step PSIM lifecycle harness | `gpu_resident_hpsi` 10.05 s at 512/1 | - | - | Adds an `NSTEPS=2` harness plus per-case copy-row extraction; all cases finish with the same two-step energy, and the profile confirms the next boundary is still `PSI0`/`HPSI`/ADDPRO-style residency rather than another immediate PSIM-only promotion. |
| `opsi-present-consumers-20260601-*` | OPSI present-or-copy consumer cleanup | `gpu_resident_hpsi_opsi` 10.13 s at 512/1 | - | - | Converts downstream projection/ADDPRO/ADDOPSI data regions to `present_or_copy*`; correctness is preserved, but the superwave OPSI build remains the real copy target. |
| `opsi-superwave-build-residency-20260601-final-v2` | Superwave OPSI build residency | `gpu_resident_hpsi_opsi` 10.36 s at 512/1 | - | - | Keeps superwave OPSI resident through build/mass scaling with one post-mass host snapshot; energy-valid and removes ADDPRO-OPSI in/out copies. |
| `superwave-projection-residency-20260601-512-nstep2-*` | Superwave projection residency | `gpu_resident_hpsi_opsi` 10.61 s at 512/1 | - | `gpu_resident_hpsi_opsi` 15.14 s at 512/4 | Enables the resident cuBLAS projection path for superwave OPSI with gamma correction; removes the post-mass OPSI host snapshot. |
| `psi0-hpsi-copy-boundaries-20260601-512-nstep2-*` | HPSI/PSI0 ETOT residency | `gpu_resident` 11.20 s at 512/1; `gpu_resident_hpsi_opsi` lowest copy | - | `gpu_resident_hpsi_opsi` 15.33 s at 512/4 | Keeps `PSI0` present from HPSI into the immediate expectation/Hamiltonian overlaps; energy-valid and removes one more 0.1345 GB copy block from HPSI/OPSI diagnostics. |
| `force-to-hpsi-psi0-residency-20260601-512-nstep2-*` | Force-to-HPSI `PSI0` residency | `gpu_resident_hpsi_opsi` 11.03 s at 512/1 | - | `gpu_resident_hpsi_opsi` 15.06 s at 512/4 | Reuses the force-loop `PSI0` device copy in the following HPSI path; energy-valid and removes the HPSI-side `PSI0` copy. |
| `vpsi-device-finish-residency-20260601-512-nstep2-*` | VPSI HPSI device finish | `gpu_resident_hpsi_opsi` 11.11 s at 512/1 | - | `gpu_resident_hpsi` 15.16 s at 512/4 | Moves the still-required HPSI transfer from the ADDPRO consumer boundary to the VPSI producer boundary; energy-valid and keeps HPSI present for ADDPRO. |
| `si64_bands-nvhpc-refresh-20260601-85882ef-1024-nstep1` | Current full matrix after HPSI/VPSI residency | `gpu_resident_hpsi_opsi` 12.83 s | `cpu` 73.32 s, `nvhpc_cpu` 69.50 s | `cpu` 167.49 s, `nvhpc_cpu` 166.88 s | Confirms wavefunction residency dominates on Spark; all-library paths remain diagnostic-only. |
| `hpsi-opsi-combo-cases-20260601-512-nstep2-*` | Combined HPSI/OPSI diagnostic keywords | `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` 9.58 s at 512/1 | - | `gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum` 15.02 s at 512/4 | Adds harness cases for HPSI+OPSI with PROJ/off-site/DENMAT combinations; all cases are energy-valid, so future standard sweeps can compare the full stack directly. |
| `si64_bands-nvhpc-standard-20260601-4fbe2cd-1024-nstep1` | Current standard refresh with combined HPSI/OPSI cases | `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` 12.21 s | `cpu` 77.79 s, `nvhpc_cpu` 75.23 s | `cpu` 672.52 s, `nvhpc_cpu` 392.07 s | Full focused residency stack is now the best 1024/1 case; eight-rank CPU is a poor resource comparison for this small smoke. |
| `si64_bands-focus-20260601-0a43a65-2048-nstep1-1r` | 2048-band focused stack validation | `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` 35.28 s | - | - | Confirms the focused stack also wins at 2048/1, so expose it as a short benchmark keyword. |
| `vpsi-boundary-20260601-d3cd6cc-512-nstep1` / `vpsi-boundary-20260601-7ad2625-smoke` | VPSI/HPSI boundary harness | `gpu_resident_stack` 6.56 s at 512/1, 28.46 s at 512/4 | - | - | Adds a focused producer-boundary harness plus seconds-sorted `PAW_VPSI_*` rows; the final smoke shows VPSI time is almost entirely GTOR/RTOG. |
| `psi0-prinfo-spark-20260601-074852` / `psi0-prinfo-terok-20260601-074852` | PSI0-to-PRINFO residency | `gpu_resident_stack` 12.11 s on Spark, 12.77 s on Terok | - | - | Keeps `PSI0` resident through `PRINFO/WRITEPDOS`, removing one more 0.1210 GB wavefunction copy; energy-valid on both systems, wall time neutral/noisy. |
| `hpsi-prop-spark-20260601-100735` / `hpsi-prop-terok-20260601-101036` | HPSI-to-propagate residency diagnostic | `gpu_resident_stack_hpsi_prop_psim_phase` 12.13 s on Spark, 12.59 s on Terok | - | - | Keeps `HPSI` resident from ETOT overlap/Hamiltonian into GPU PSIM propagation; removes the extra 0.1210 GB `PROP_HPSI_IN` copy from the PSIM-phase diagnostic and stays energy-valid on both systems. |
| `psim-switch-accdelete-spark-20260601-104028` / `psim-switch-accdelete-terok-20260601-104027` | NSTEPS=2 PSIM switch residency | `gpu_resident_stack_hpsi_prop_psim_switch` 20.61 s on Spark, 20.25 s on Terok | - | - | Carries resident `PSIM` through `WAVES$SWITCH` as the next `PSI0`, removes one more 0.1210 GB force-side `PSI0` copy per two-step run, and fixes the Spark partially-present delete failure. |
| `accdims-switch-spark-20260601-105947` / `accdims-switch-terok-20260601-105946` | Resident wavefunction dimension tracking | `gpu_resident_stack_hpsi_prop_psim_switch` 20.84 s on Spark, 23.50 s on Terok | - | 4-rank smokes OK | Stores resident `PSI0`/`PSIM`/`HPSI` dimensions in `WVSET_TYPE` and routes lifecycle cleanup through mark/clear/delete helpers, preserving the previous switch-residency transfer pattern. |
| `vpsi-cufft-refresh-spark-20260601-110633` / `vpsi-cufft-refresh-terok-20260601-110633` | Current VPSI/cuFFT refresh | `gpu_resident_stack_cufft` neutral/slightly favorable; force cuFFT slower | - | - | Rechecks cuFFT after the latest residency work. Threshold-gated cuFFT remains harmless, but forced cuFFT inflates transfer to 10.97/19.62 GB and slows VPSI strongly. |
| `lazy-scratch-*-20260601-1118/1121` | Lazy host scratch allocation for resident `PRO` cache paths | `gpu_resident_stack` 12.12 s on Spark, 11.99 s on Terok at 1024/1 | - | 4-rank smokes OK | Avoids building unused host `GVEC`/`PRO`/`EIGR` scratch in cached resident `PRO` projection and addproduct paths; both cache and host-PRO ablation paths stay energy-valid. |
| `psim-stack-default-*-20260601-1140/1145/1150` | Pre-lifecycle-fix PSIM stack-default probe | Serial opt-in saves 0.1204 GB at 1024/2; Terok 4-rank was much slower before the dimension/lifecycle cleanup | - | Terok 4-rank regression in old probe | Superseded by the later resident-dimension cleanup and default retest below; kept as cautionary history. |
| existing `lazy-scratch-1024-*` profiles re-summarized | FFT/VPSI benchmark collector fields | Spark `gpu_resident_stack`: `vpsi_s=1.2877`, `vpsi_gtor_s=0.6415`, `vpsi_rtog_s=0.6279` | - | tooling OK | Adds `pw_fft_gtor_s`, `pw_fft_rtog_s`, `vpsi_s`, `vpsi_gtor_s`, and `vpsi_rtog_s` to `benchmark_summary.py`; `run_vpsi_boundary.sh` now includes `PW_FFT_*` rows in its FFT-phase table. |
| `stack-default-psim-*-20260601-1206/1208` | PSIM switch promoted into stack default after lifecycle cleanup | `gpu_resident_stack` now matches the explicit switch case at 512/2: 9.73 s on Spark, 11.07 s on Terok; `NSTEPS=1` smokes OK | - | 4-rank smokes OK | `CPPAW_GPU_RESIDENCY_STACK` now enables PSIM propagation/phase/switch and HPSI-to-propagate residency by default, removing `ACC_COPY_ORTHO_PSIM_IN` and saving 0.0666 GB at 512/2. |
| `current-stack-spark-1024-nstep1-*` | Fresh Spark stack-default refresh plus CPU-build guard | `gpu_resident_stack` 12.19 s; threshold-gated cuFFT 12.19 s; forced cuFFT 15.30 s | `nvhpc_cpu` 72.41 s | `nvhpc_cpu` 166.04 s | Fixes the non-CUBLAS `WAVES$HPSI` CPU build guard and confirms the current stack default is energy-valid and much faster than the CPU references. Forced cuFFT remains diagnostic-only because it raises transfer volume to 10.05 GB. |
| `dual-switch-*-20260601-1249/1250/1252` | Bidirectional PSIM/PSI0 switch residency | `gpu_resident_stack` keeps the old `PSI0` as resident `PSIM` across `WAVES$SWITCH`; `ACC_COPY_PROP_PSIM_IN` drops from 3 calls to 1 at 1024/3 | - | 4-rank smoke OK | Removes 0.2421 GB of repeated propagation input traffic in the 1024/3 smoke, while the 1024/1 and 512/1x4 smokes remain energy-valid. |
| `serial3dfft-accmap-final-*-20260601-*` | Device-side sparse/full-grid mapping for single-rank 3D cuFFT | Terok/A40: `gpu_resident_stack_serial3dfft_accmap` 10.62 s at NSTEPS=1 and 25.55 s at NSTEPS=3 | Spark/GB10: regular `gpu_resident_stack_serial3dfft` 11.02 s at NSTEPS=1 and 23.44 s at NSTEPS=3 | - | Adds an explicit diagnostic case for device-side mapping. It helps on A40 but hurts on GB10, so it remains opt-in via `CPPAW_FFT_SERIAL_3D_ACC_MAP=1`. |
| `setup-psim-isolated-*-20260601-*` | Setup `PSIM` residency in the focused stack | Spark/GB10: `gpu_resident_stack` 27.96 s at NSTEPS=3; Terok/A40: 29.54 s at NSTEPS=3 | Ablation with `CPPAW_GPU_SETUP_PSIM_RESIDENCY=0`: 28.23 s on Spark and 39.17 s on Terok at NSTEPS=3 | NSTEPS=1 energy-valid on both systems | Keeps initial setup `PSIM` resident into Gram-Schmidt/propagation; removes one 0.1210 GB propagation copy without broadening non-phase PSIM diagnostics. |
| `accmap-present-nstep*-spark/terok` | Present-input reuse in the serial 3D cuFFT ACCMAP path | Transfer estimate drops from 4.0171 to 3.6540 GB at NSTEPS=1 and from 10.5201 to 9.6729 GB at NSTEPS=3 | Wall time remains mixed: Terok NSTEPS=3 improves to 24.71 s, Spark NSTEPS=3 is 41.99 s | NSTEPS=1 energy-valid on both systems | Changes the ACCMAP data region to `PRESENT_OR_COPYIN` for the FFT input vector. This reduces redundant host-to-device traffic but does not make ACCMAP a Spark default. |
| `hpsi-rtog-*-20260601-*` | HPSI RTOG output residency in the serial 3D FFT ACCMAP path | Spark/GB10: output-present accounting lowers the validated NSTEPS=1 transfer estimate to 3.4119 GB; wall time is modestly favorable in the initial sweep | Terok/A40: NSTEPS=1 remains mixed/noisy, but the transfer accounting is identical | NSTEPS=1 energy-valid on both systems | Adds opt-in `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1` / `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog`; useful as a producer-boundary diagnostic, not a default promotion. |
| `vpsi-internal-*-20260601-*` | Resident `WAVES_VPSI` real-space scratch in the serial 3D FFT ACCMAP path | Spark/GB10: isolated scratch residency lowers transfer to 2.2890 GB at NSTEPS=1; combined with HPSI-RTOG it reaches 2.0469 GB | Terok/A40: isolated scratch residency lowers transfer identically, but combined HPSI-RTOG is the better diagnostic at NSTEPS=3 | NSTEPS=1 energy-valid on both systems | Adds separate `gpu_resident_stack_serial3dfft_accmap_vpsi_internal` and combined `*_hpsi_rtog_vpsi_internal` cases; useful and measurable, but still opt-in because wall time is case/system noisy. |
| `accmap-cache-*-20260601-*` | Cached serial 3D ACCMAP work/map arrays | Spark/GB10: combined HPSI-RTOG+VPSI internal cache improves NSTEPS=3 from 39.09 to 34.58 s and lowers transfer from 4.8517 to 4.5655 GB | Terok/A40: same cache improves NSTEPS=3 from 33.10 to 23.66 s and lowers transfer identically | NSTEPS=1 energy-valid on both systems | Adds opt-in `CPPAW_FFT_SERIAL_3D_ACC_CACHE=1` / `*_hpsi_rtog_vpsi_internal_cache`; validated as a useful diagnostic, but still not a default because the serial ACCMAP path remains system-dependent. |
| `accmap-cleanup-final-*-20260601-*` | ACCMAP cache cleanup and non-CUBLAS build guard | Spark/GB10: final cache smoke is 8.16 s at 512/1 and 6.53 s at 256/4 | Terok/A40: final cache smoke is 5.25 s at 512/1 and 8.58 s at 256/4 | Energy-valid; `nvhpc_profile` and GPU serial/parallel builds pass on both systems | Moves the cached ACCMAP state into `PLANEWAVE_MODULE`, releases it through `PLANEWAVE$ACC_CLEANUP`, and restores the non-CUBLAS `nvhpc_profile` build by guarding setup-PSIM residency code. |
| `auto-standard-spark-20260601-164010` / `auto-standard-terok-fftw-20260601-164730` | Capability-driven standard smoke | Spark/GB10: `gpu_resident_stack_serial3dfft` 11.13 s, `gpu_resident_stack` 12.14 s | Spark: `cpu` 73.51 s, `nvhpc_cpu` 69.11 s; 8-rank CPU 166.30/166.41 s | Terok x86_64 builds and runs after local FFTW plus NVHPC compiler BLAS/LAPACK fallback | Confirms the auto recommendation path and the new host-library gating. Spark remains the performance reference; Terok is the x86/NVHPC portability check. |
| `si64-bands2048-focused-spark-20260601-174042` / `si64-bands2048-focused-terok-20260601-174041` | 2048-band resource comparison | Spark: `gpu_resident_stack_serial3dfft_force_dedpro` 32.06 s; Terok: `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache` 32.79 s | Spark `nvhpc_cpu` 388.98 s; Terok `nvhpc_cpu` 1109.15 s | Spark `nvhpc_cpu` 1122.63 s; Terok `nvhpc_cpu` 1341.45 s | Larger band stress confirms that one MPI rank plus one GPU beats both one-rank and eight-rank CPU/NVHPC decisively; ACCMAP remains system-dependent. |

The latest auto-standard smoke runs live at:

```
Spark: /home/kuehne88/cp-paw-nvhpc-auto-20260601-163927/tests/profile/si64/runs/auto-standard-spark-20260601-164010
Terok: /home/kuehne88/cp-paw-nvhpc-auto-20260601-163926/tests/profile/si64/runs/auto-standard-terok-fftw-20260601-164730
```

## 2026-06-01 Auto Standard And Host Library Gating

The standard wrapper now consumes the capability helper's recommended case
lists by default and can auto-build the required profile targets. The helper no
longer recommends GPU or NVHPC CPU cases unless the host numerical stack is
usable: host FFTW, including `fftw3.f03`, and host BLAS/LAPACK are reported
explicitly. If either is missing, `run_nvhpc_standard.sh` resolves both auto
case lists to empty, logs `SKIP all suites empty case lists`, and exits cleanly
instead of failing later inside `paw_build.sh`.

Spark C86C has the full NVHPC/GB10 stack available from the system SDK. The
auto-standard run built all required targets and completed the GPU, one-rank
CPU, and eight-rank CPU suites:

| Case | Ranks | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident_stack_serial3dfft` | 1 | 11.13 s | 5.4859 GB | yes |
| `gpu_resident_stack` | 1 | 12.14 s | 1.3892 GB | yes |
| `gpu_resident_stack_cufft` | 1 | 12.47 s | 1.3892 GB | yes |
| `gpu_resident_off` | 1 | 43.77 s | 0.0000 GB | yes |
| `nvhpc_cpu` | 1 | 69.11 s | 0.0000 GB | yes |
| `cpu` | 1 | 73.51 s | 0.0000 GB | yes |
| `cpu` | 8 | 166.30 s | 0.0000 GB | yes |
| `nvhpc_cpu` | 8 | 166.41 s | 0.0000 GB | yes |

The performance conclusion is unchanged but now reproduced through the auto
path: one GPU with the resident stack is much faster than the 1-rank CPU
reference for this Si64 smoke, and also faster than the 8-rank CPU resource
comparison. `gpu_resident_stack_serial3dfft` is the fastest one-step result, but
it still carries the explicit full-grid copy estimate, so the routine default
stays the conservative resident stack while serial 3D FFT remains an opt-in
diagnostic.

Terok initially exposed the portability gap: NVHPC 24.5 on x86_64 provides
cuBLAS, cuFFT, cuSOLVER, cuTENSOR, NCCL, and NVSHMEM, but no usable host NVPL
FFTW library. Installing FFTW 3.3.10 into
`/home/kuehne88/opt/fftw-3.3.10` and exporting its `PKG_CONFIG_PATH` lets the
helper report:

```
host_fftw=yes path=pkg-config:fftw3 prefix=/home/kuehne88/opt/fftw-3.3.10
host_blas_lapack=yes path=/home/kuehne88/opt/nvidia/hpc_sdk/Linux_x86_64/2024/compilers/lib/lib{blas,lapack}.so
```

The corresponding build fallback adds the NVHPC compiler `libblas.so` and
`liblapack.so` when NVPL is absent. With that setup the same auto-standard
matrix builds and runs on Terok:

| Case | Ranks | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | --- |
| `gpu_resident_stack_serial3dfft` | 1 | 11.95 s | 5.4859 GB | yes |
| `gpu_resident_stack_cufft` | 1 | 12.17 s | 1.3892 GB | yes |
| `gpu_resident_stack` | 1 | 12.71 s | 1.3892 GB | yes |
| `gpu_resident_off` | 1 | 77.75 s | 0.0000 GB | yes |
| `nvhpc_cpu` | 1 | 190.96 s | 0.0000 GB | yes |
| `cpu` | 1 | 193.34 s | 0.0000 GB | yes |
| `cpu` | 8 | 142.11 s | 0.0000 GB | yes |
| `nvhpc_cpu` | 8 | 142.41 s | 0.0000 GB | yes |

The Terok timings were taken while other CP2K/GauXC work was visible on the
machine, so they are a portability and correctness check rather than a clean
performance comparison. Still, every selected case is energy-valid, and the
same capability-driven case matrix now works on both aarch64/GB10 with NVHPC
26.3 and x86_64/A40 with NVHPC 24.5 plus a host FFTW install.

## 2026-06-01 HPSI RTOG Output Residency Diagnostic

The serial 3-D FFT ACCMAP follow-up can create `HPSI` on the GPU before
`WAVES_VPSI`, let the RTOG mapping write into that device copy, and then skip
the old `ACC_COPY_VPSI_HPSI_IN` / update boundary when the FFT path actually
used ACCMAP. It is gated by `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1` and exposed
as `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog`.

Run directories:

```
runs/hpsi-rtog-spark-rep3-20260601-150822-{base,rtog}
runs/hpsi-rtog-terok-rep3-20260601-150821-{base,rtog}
runs/hpsi-rtog-spark-n3-20260601-151020-{base,rtog}
runs/hpsi-rtog-terok-n3-20260601-151019-{base,rtog}
runs/hpsi-rtog-accounting-spark-20260601-151835
runs/hpsi-rtog-accounting-terok-20260601-151834
```

| System | NSTEPS | Case | Wall time | VPSI time | Transfer estimate | Energy check |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| Spark GB10 | 1 | ACCMAP baseline, 3-run avg | 16.37 s | 1.2115 s | 3.6540 GB | yes |
| Spark GB10 | 1 | `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1`, 3-run avg | 15.61 s | 0.2467 s | 3.5330 GB | yes |
| Spark GB10 | 3 | ACCMAP baseline | 41.60 s | 3.7410 s | 9.6729 GB | n/a |
| Spark GB10 | 3 | `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1` | 38.72 s | 0.7037 s | 9.3098 GB | n/a |
| Terok A40 | 1 | ACCMAP baseline, 3-run avg | 10.52 s | 0.5596 s | 3.6540 GB | yes |
| Terok A40 | 1 | `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1`, 3-run avg | 13.67 s | 0.4329 s | 3.5330 GB | yes |
| Terok A40 | 3 | ACCMAP baseline, noisy single run | 43.16 s | 9.4256 s | 9.6729 GB | n/a |
| Terok A40 | 3 | `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1`, noisy single run | 25.12 s | 1.8182 s | 9.3098 GB | n/a |

The profile row changes as intended: the baseline reports
`ACC_COPY_VPSI_HPSI_IN` with 0.1210 GB at NSTEPS=1, while the opt-in path
reports `ACC_PRESENT_VPSI_HPSI` and removes that transfer. Spark benefits in
both short smokes. Terok is noisy and the NSTEPS=1 average is slower, so this
stays an explicit diagnostic for broader producer/consumer residency work rather
than part of `CPPAW_GPU_RESIDENCY_STACK`.

A follow-up accounting-only correction also subtracts ACCMAP outputs that are
already present on the device. The earlier HPSI-RTOG transfer estimates above
therefore overstate the current instrumented path by one resident `HPSI` output
vector, 0.1210 GB per step. Current NSTEPS=1 spot checks give:

| System | Case | Wall time | VPSI time | Transfer estimate | `ACC_COPY_SERIAL3D_ACC_MAP` | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | HPSI-RTOG output-present accounting | 15.83 s | 0.2495 s | 3.4119 GB | 2.1438 GB | yes |
| Terok A40 | HPSI-RTOG output-present accounting | 11.36 s | 1.0396 s | 3.4119 GB | 2.1438 GB | yes |

## 2026-06-01 VPSI Internal Scratch Residency Diagnostic

The next producer-internal diagnostic keeps the scalar-spin `WAVES_VPSI`
real-space scratch `PSIOFR` present on the GPU between GTOR, the local potential
multiplication, and RTOG. It is deliberately narrow: `NDIM=1`, HPSI residency
active, serial 3-D FFT ACCMAP active, and `CPPAW_GPU_VPSI_INTERNAL_RESIDENCY=1`.
The benchmark case is
`gpu_resident_stack_serial3dfft_accmap_vpsi_internal`; the combined
`gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal` case also
enables `CPPAW_GPU_VPSI_HPSI_RTOG_RESIDENCY=1`.

Run directories:

```
Spark: tests/profile/si64/runs/vpsi-internal-spark-20260601-152856
Spark: tests/profile/si64/runs/vpsi-internal-n3-spark-20260601-153002
Spark: tests/profile/si64/runs/vpsi-internal-isolated-spark-20260601-153725
Spark: tests/profile/si64/runs/vpsi-internal-isolated-n3-spark-20260601-153837
Terok: tests/profile/si64/runs/vpsi-internal-terok-20260601-152855
Terok: tests/profile/si64/runs/vpsi-internal-n3-terok-20260601-153002
Terok: tests/profile/si64/runs/vpsi-internal-isolated-terok-20260601-153724
Terok: tests/profile/si64/runs/vpsi-internal-isolated-n3-terok-20260601-153836
```

| System | NSTEPS | Case | Wall time | VPSI time | Transfer estimate | Energy check |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| Spark GB10 | 1 | HPSI-RTOG ACCMAP | 16.96 s | 0.2461 s | 3.4119 GB | yes |
| Spark GB10 | 1 | HPSI-RTOG + VPSI internal scratch | 15.49 s | 0.1057 s | 2.0469 GB | yes |
| Spark GB10 | 3 | HPSI-RTOG ACCMAP | 43.19 s | 0.7854 s | 8.9467 GB | n/a |
| Spark GB10 | 3 | HPSI-RTOG + VPSI internal scratch | 38.10 s | 0.3174 s | 4.8517 GB | n/a |
| Terok A40 | 1 | HPSI-RTOG ACCMAP | 10.81 s | 0.6311 s | 3.4119 GB | yes |
| Terok A40 | 1 | HPSI-RTOG + VPSI internal scratch | 10.24 s | 0.1327 s | 2.0469 GB | yes |
| Terok A40 | 3 | HPSI-RTOG ACCMAP | 32.70 s | 5.0492 s | 8.9467 GB | n/a |
| Terok A40 | 3 | HPSI-RTOG + VPSI internal scratch | 23.82 s | 0.4196 s | 4.8517 GB | n/a |

Representative profile rows are identical on Spark and Terok for the one-step
case: `ACC_COPY_SERIAL3D_ACC_MAP` drops from 2.1438 GB to 0.7782 GB,
`ACC_PRESENT_VPSI_PSIOFR` records 576 resident scratch uses, and
`ACC_COPY_VPSI_V_IN` adds only 0.0006 GB for the local potential input. This is
therefore a real transfer reduction, not just a row-label change.

The isolation follow-up splits the effects: `CPPAW_GPU_VPSI_INTERNAL_RESIDENCY`
alone keeps `PSIOFR` resident but still leaves the old `ACC_COPY_VPSI_HPSI_IN`
boundary, while the combined case also makes RTOG produce resident `HPSI`:

| System | NSTEPS | Case | Wall time | VPSI time | Transfer estimate | Energy check |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| Spark GB10 | 1 | ACCMAP baseline | 16.94 s | 1.2259 s | 3.6540 GB | yes |
| Spark GB10 | 1 | VPSI internal scratch only | 15.06 s | 0.8582 s | 2.2890 GB | yes |
| Spark GB10 | 1 | VPSI internal scratch + HPSI-RTOG | 17.42 s | 0.1062 s | 2.0469 GB | yes |
| Spark GB10 | 3 | ACCMAP baseline | 38.17 s | 3.2642 s | 9.6729 GB | n/a |
| Spark GB10 | 3 | VPSI internal scratch only | 37.27 s | 2.6285 s | 5.5779 GB | n/a |
| Spark GB10 | 3 | VPSI internal scratch + HPSI-RTOG | 35.02 s | 0.3156 s | 4.8517 GB | n/a |
| Terok A40 | 1 | ACCMAP baseline | 11.08 s | 0.7460 s | 3.6540 GB | yes |
| Terok A40 | 1 | VPSI internal scratch only | 10.51 s | 0.2843 s | 2.2890 GB | yes |
| Terok A40 | 1 | VPSI internal scratch + HPSI-RTOG | 10.38 s | 0.1321 s | 2.0469 GB | yes |
| Terok A40 | 3 | ACCMAP baseline | 24.77 s | 2.0431 s | 9.6729 GB | n/a |
| Terok A40 | 3 | VPSI internal scratch only | 37.32 s | 1.0991 s | 5.5779 GB | n/a |
| Terok A40 | 3 | VPSI internal scratch + HPSI-RTOG | 24.47 s | 0.4335 s | 4.8517 GB | n/a |

The isolated case therefore proves the `PSIOFR` residency accounting and
correctness independently, but the combined case is the more useful ACCMAP
diagnostic. Neither should be folded into `CPPAW_GPU_RESIDENCY_STACK` yet,
because the serial 3-D ACCMAP path itself is still system-dependent.

## 2026-06-01 Serial 3D ACCMAP Phase Profiling

The ACCMAP path now reports the visible device-side mapping subphases separately:
`ACC_SERIAL3D_GTOR_ZERO`, `ACC_SERIAL3D_GTOR_SCATTER`,
`ACC_SERIAL3D_GTOR_GATHER`, `ACC_SERIAL3D_RTOG_LOAD`, and
`ACC_SERIAL3D_RTOG_GATHER`. The cuFFT call itself is already reported by
`CUFFT3D_C8_PRESENT`.

Run directories:

```
Spark: tests/profile/si64/runs/accmap-phase-profile-spark-20260601-154539
Terok: tests/profile/si64/runs/accmap-phase-profile-terok-20260601-154539
```

Both systems used `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal`
with `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`, `RANKS=1`.

| System | Wall time | FFT time | VPSI time | Mapping kernels | `CUFFT3D_C8_PRESENT` | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | 15.48 s | 5.1413 s | 0.1069 s | 0.0919 s | 0.0393 s | 2.0469 GB | yes |
| Terok A40 | 10.03 s | 0.5504 s | 0.1359 s | 0.0970 s | 0.0431 s | 2.0469 GB | yes |

The visible mapping kernels and cuFFT call are much smaller than the reported
FFT envelope, especially on Spark. This points away from the map kernels
themselves and toward per-call OpenACC data-region/runtime overhead and repeated
temporary full-grid workspace lifetimes as the next serial-3D ACCMAP target.

## 2026-06-01 Serial 3D ACCMAP Work/Map Cache

The immediate follow-up keeps the serial 3-D ACCMAP full-grid `WORK` array and
map arrays present across calls when `CPPAW_FFT_SERIAL_3D_ACC_CACHE=1` is set.
The benchmark case
`gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache` combines
that cache with HPSI RTOG output residency and resident `WAVES_VPSI` real-space
scratch. This is intentionally opt-in because serial ACCMAP remains
system-dependent, but the cached device allocations now have an explicit
plane-wave accelerator cleanup hook. The case tests whether repeated temporary
workspace creation and map transfers explain the remaining ACCMAP overhead.

Run directories:

```
Spark: tests/profile/si64/runs/accmap-cache-spark-20260601-155300
Spark: tests/profile/si64/runs/accmap-cache-n3-spark-20260601-155518
Terok: tests/profile/si64/runs/accmap-cache-terok-20260601-155301
Terok: tests/profile/si64/runs/accmap-cache-n3-terok-20260601-155616
Spark 4-rank smoke: tests/profile/si64/runs/accmap-cache-parallel-smoke-spark-20260601-160305
Terok 4-rank smoke: tests/profile/si64/runs/accmap-cache-parallel-smoke-terok-20260601-160304
```

| System | NSTEPS | Case | Wall time | FFT time | VPSI time | Transfer estimate | `ACC_COPY_SERIAL3D_ACC_MAP` | Energy check |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | 1 | HPSI-RTOG + VPSI internal | 15.46 s | 5.3818 s | 0.1071 s | 2.0469 GB | 0.7782 GB | yes |
| Spark GB10 | 1 | HPSI-RTOG + VPSI internal + cache | 15.20 s | 5.3771 s | 0.0954 s | 1.9516 GB | 0.6828 GB | yes |
| Spark GB10 | 3 | HPSI-RTOG + VPSI internal | 39.09 s | 17.2552 s | 0.3199 s | 4.8517 GB | 2.3347 GB | n/a |
| Spark GB10 | 3 | HPSI-RTOG + VPSI internal + cache | 34.58 s | 13.0891 s | 0.2859 s | 4.5655 GB | 2.0484 GB | n/a |
| Terok A40 | 1 | HPSI-RTOG + VPSI internal | 10.67 s | 0.6057 s | 0.1350 s | 2.0469 GB | 0.7782 GB | yes |
| Terok A40 | 1 | HPSI-RTOG + VPSI internal + cache | 10.68 s | 0.7889 s | 0.1189 s | 1.9516 GB | 0.6828 GB | yes |
| Terok A40 | 3 | HPSI-RTOG + VPSI internal | 33.10 s | 5.6455 s | 3.0183 s | 4.8517 GB | 2.3347 GB | n/a |
| Terok A40 | 3 | HPSI-RTOG + VPSI internal + cache | 23.66 s | 1.6763 s | 0.3509 s | 4.5655 GB | 2.0484 GB | n/a |

The cache adds one `ACC_COPY_SERIAL3D_ACC_MAP_CACHE` row of only 0.000055 GB
for this Si64 grid and one `ACC_CREATE_SERIAL3D_ACC_WORK` row, then reduces the
repeated map-copy estimate by 0.0954 GB at NSTEPS=1 and 0.2863 GB at NSTEPS=3.
Spark sees a small one-step gain and a clear three-step gain. Terok is neutral
at one step but strongly favorable at three steps in this run. The result is
therefore strong enough to keep the cache as a validated diagnostic switch, but
not strong enough to make serial 3-D ACCMAP itself part of the recommended
default stack.

## 2026-06-01 Serial 3D ACCMAP Cache Cleanup

The cache follow-up moves the saved ACCMAP work array and cache bookkeeping into
`PLANEWAVE_MODULE` and releases the OpenACC data in `PLANEWAVE$ACC_CLEANUP`,
which is called after the timing/profiling report. Non-`CPPVAR_CUFFT_ACC` builds
compile the cleanup routine as a no-op. The same patch also guards the setup
`PSIM` residency block in `WAVES$GRAMMSCHMIDT`, restoring the plain
`nvhpc_profile` build without `CPPVAR_CUBLAS_ACC`.

Build validation:

```
Spark: nvhpc_profile, nvhpc_gpu_acc_residency_profile, nvhpc_gpu_acc_residency_profile_parallel
Terok: nvhpc_profile, nvhpc_gpu_acc_residency_profile, nvhpc_gpu_acc_residency_profile_parallel
```

Runtime smoke directories:

```
Spark serial: tests/profile/si64/runs/accmap-cleanup-final-serial-spark-20260601-161939
Spark 4-rank: tests/profile/si64/runs/accmap-cleanup-final-4rank-spark-20260601-161947
Terok serial: tests/profile/si64/runs/accmap-cleanup-final-serial-terok-20260601-161938
Terok 4-rank: tests/profile/si64/runs/accmap-cleanup-final-4rank-terok-20260601-161944
```

| System | Ranks | Empty bands | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | 1 | 512 | 8.16 s | 1.0923 GB | yes |
| Spark GB10 | 4 | 256 | 6.53 s | 0.6901 GB | yes |
| Terok A40 | 1 | 512 | 5.25 s | 1.0923 GB | yes |
| Terok A40 | 4 | 256 | 8.58 s | 0.6901 GB | yes |

Interpretation: this is a robustness patch for the diagnostic ACCMAP cache, not
a default-promotion signal. The cache now has a real cleanup boundary and the
CPU-only NVHPC profile build remains valid.

## 2026-06-01 Current Stack Default Refresh

After promoting PSIM switch/phase/propagation residency into
`CPPAW_GPU_RESIDENCY_STACK`, Spark C86C was retested from a fresh checkout using
`TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`, and three repeats per case.
The same validation also caught and fixed a CPU-build regression where
`WAVES$HPSI` called the HPSI residency-clear helper outside the
`CPPVAR_CUBLAS_ACC` guard.

Run directories:

```
runs/current-stack-spark-1024-nstep1-20260601-121803
runs/current-stack-spark-1024-nstep1-cpu1-20260601-122316
runs/current-stack-spark-1024-nstep1-cpu8-20260601-122700
```

| Suite | Case | Ranks | Repeats | Wall time | Transfer estimate | Energy check | Interpretation |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| GPU | `gpu_resident_stack` | 1 | 3 | 12.19 s | 1.51 GB | yes | Current recommended stack keyword; broad residency remains the decisive lever. |
| GPU | `gpu_resident_stack_cufft` | 1 | 3 | 12.19 s | 1.51 GB | yes | Threshold-gated cuFFT is neutral/harmless for this case. |
| GPU | `gpu_resident_stack_cufft_force` | 1 | 3 | 15.30 s | 10.05 GB | yes | Forced cuFFT is slower and increases transfer volume strongly; keep diagnostic-only. |
| CPU | `nvhpc_cpu` | 1 | 3 | 72.41 s | 0.00 GB | yes | Best direct one-rank CPU reference in this refresh. |
| CPU | `nvhpc_cpu` | 8 | 3 | 166.04 s | 0.00 GB | yes | Poor resource comparison for this small smoke; MPI/setup overhead dominates. |

The main remaining copy rows in `gpu_resident_stack` are now the ZGEMM matrix
copy-out, setup/Gram PSIM boundaries, `VPSI` HPSI input, propagation PSIM input,
and the final orthogonalization PSIM output. That reinforces the current design
direction: keep widening wavefunction/projector residency across producer and
consumer boundaries; do not promote forced cuFFT for this workload.

## 2026-06-01 Bidirectional PSIM/PSI0 Switch Residency

The previous stack default carried resident `PSIM` through `WAVES$SWITCH` as the
next `PSI0`, but it deleted the old resident `PSI0`. After the pointer swap that
old `PSI0` is exactly the next-step `PSIM`, so the following propagation copied
it back from the host. The switch now preserves both resident sides when
`CPPAW_GPU_PSIM_SWITCH_RESIDENCY` is enabled: old `PSIM` becomes resident
`PSI0`, and old `PSI0` becomes resident `PSIM`.

Run directories:

```
runs/psim-copy-scaling-20260601-124358
runs/energy-triage-20260601-124513
runs/dual-switch-residency-20260601-124939
runs/dual-switch-smoke-20260601-125043
runs/dual-switch-4rank-smoke-20260601-125215
```

| Case | Build | Ranks | NSTEPS | Wall time | Transfer estimate | `ACC_COPY_PROP_PSIM_IN` | Energy check | Interpretation |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Before patch | serial GPU stack | 1 | 3 | 28.08 s | 3.2414 GB | 3 calls, 0.3631 GB | disabled | Baseline repeatedly copied the previous-step `PSIM` into propagation. |
| After patch | serial GPU stack | 1 | 3 | 27.97 s | 2.9994 GB | 1 call, 0.1210 GB | disabled | Only the first propagation needs the host `PSIM`; two later steps record `ACC_PRESENT_PROP_PSIM`. |
| After patch | serial GPU stack | 1 | 1 | 12.10 s | 1.5102 GB | 1 call, 0.1210 GB | yes | Standard 1024-band smoke remains energy-valid. |
| After patch | parallel GPU stack | 4 | 1 | 9.81 s | 1.2401 GB | 4 calls, 0.0672 GB | yes | Parallel build and 4-rank smoke remain energy-valid. |

The `NSTEPS=3` energy value is intentionally not checked against the one-step
Si64 reference: a GPU-off run in the same patched binary gives the identical
final value (`208.886424`), so this is benchmark-harness behavior for the
multi-step wavefunction dynamics, not a GPU correctness regression. The harness
now applies the built-in Si64 energy reference only to `NSTEPS=1` unless an
explicit `EXPECTED_ENERGY` is supplied.

## 2026-06-01 Serial 3D cuFFT Device Mapping

The first single-rank 3D cuFFT diagnostic still copied the full temporary
`WORK(NR1,NR2,NR3)` grid through the generic `LIB$3DFFTC8` wrapper. A follow-up
diagnostic keeps the sparse/full-grid mapping on the GPU and runs cuFFT on a
present full-grid buffer. It is selected only when both
`CPPAW_FFT_SERIAL_3D=1` and `CPPAW_FFT_SERIAL_3D_ACC_MAP=1` are set; the
benchmark keyword is `gpu_resident_stack_serial3dfft_accmap`.

Run directories:

```
Spark: tests/profile/si64/runs/serial3dfft-accmap-final-spark-20260601-133950
Spark: tests/profile/si64/runs/serial3dfft-accmap-final-spark-nsteps3-20260601-134222
Terok: tests/profile/si64/runs/serial3dfft-accmap-final-terok-20260601-133949
Terok: tests/profile/si64/runs/serial3dfft-accmap-final-terok-nsteps3-20260601-134222
```

| System | Case | NSTEPS | Wall time | FFT time | VPSI time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | `gpu_resident_stack` | 1 | 12.37 s | 3.7735 s | 1.3191 s | 1.5102 GB | yes |
| Spark GB10 | `gpu_resident_stack_serial3dfft` | 1 | 11.02 s | 0.9500 s | 0.3643 s | 5.6070 GB | yes |
| Spark GB10 | `gpu_resident_stack_serial3dfft_accmap` | 1 | 18.38 s | 7.9962 s | 1.5265 s | 4.0171 GB | yes |
| Terok A40 | `gpu_resident_stack` | 1 | 12.27 s | 3.4416 s | 1.4953 s | 1.5102 GB | yes |
| Terok A40 | `gpu_resident_stack_serial3dfft` | 1 | 11.54 s | 2.7489 s | 1.0131 s | 5.6070 GB | yes |
| Terok A40 | `gpu_resident_stack_serial3dfft_accmap` | 1 | 10.62 s | 0.9957 s | 0.3850 s | 4.0171 GB | yes |
| Spark GB10 | `gpu_resident_stack` | 3 | 28.09 s | 11.2278 s | 3.9124 s | 2.9994 GB | n/a |
| Spark GB10 | `gpu_resident_stack_serial3dfft` | 3 | 23.44 s | 2.6494 s | 1.0490 s | 15.2897 GB | n/a |
| Spark GB10 | `gpu_resident_stack_serial3dfft_accmap` | 3 | 41.22 s | 19.0376 s | 3.7872 s | 10.5201 GB | n/a |
| Terok A40 | `gpu_resident_stack` | 3 | 32.86 s | 9.5358 s | 3.8142 s | 2.9994 GB | n/a |
| Terok A40 | `gpu_resident_stack_serial3dfft` | 3 | 35.39 s | 13.4791 s | 5.8015 s | 15.2897 GB | n/a |
| Terok A40 | `gpu_resident_stack_serial3dfft_accmap` | 3 | 25.55 s | 2.8725 s | 1.5215 s | 10.5201 GB | n/a |

The result is intentionally mixed. Device-side mapping reduces transfer volume
against the original serial-3D cuFFT wrapper, but Spark GB10 spends much more
time in the mapping kernels. Terok A40 benefits clearly. The path therefore
stays a separate diagnostic case rather than replacing
`gpu_resident_stack_serial3dfft`.

A later follow-up lets the ACCMAP path reuse a resident input vector instead of
forcing a fresh `COPYIN` for every `PLANEWAVE$FFT` call. The output side is still
left as `COPYOUT` for conservative host visibility when no outer device mapping
exists. On the same 1024-empty-band Si64 probe, the profile-estimated transfer
volume drops from 4.0171 GB to 3.6540 GB at `NSTEPS=1` and from 10.5201 GB to
9.6729 GB at `NSTEPS=3`. The wall-time outcome remains mixed:

| System | Case | NSTEPS | Wall time | FFT time | VPSI time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Spark GB10 | `gpu_resident_stack_serial3dfft_accmap` with present input | 1 | 16.79 s | 6.2755 s | 1.3117 s | 3.6540 GB | yes |
| Spark GB10 | `gpu_resident_stack_serial3dfft_accmap` with present input | 3 | 41.99 s | 19.3489 s | 3.7334 s | 9.6729 GB | n/a |
| Terok A40 | `gpu_resident_stack_serial3dfft_accmap` with present input | 1 | 11.49 s | 1.5243 s | 0.9400 s | 3.6540 GB | yes |
| Terok A40 | `gpu_resident_stack_serial3dfft_accmap` with present input | 3 | 24.71 s | 2.7637 s | 1.7184 s | 9.6729 GB | n/a |

Run directories:

```
Spark: tests/profile/si64/runs/accmap-present-nstep1-spark
Spark: tests/profile/si64/runs/accmap-present-nstep3-spark
Terok: tests/profile/si64/runs/accmap-present-nstep1-terok
Terok: tests/profile/si64/runs/accmap-present-nstep3-terok
```

## 2026-06-01 Setup PSIM Residency

The focused stack now has a separate setup `PSIM` residency switch. The stack
default enables it through `CPPAW_GPU_RESIDENCY_STACK=1`, while the ablation
case `gpu_resident_stack_setup_psim_host` sets
`CPPAW_GPU_SETUP_PSIM_RESIDENCY=0` and leaves the rest of the stack unchanged.
The implementation is additionally gated by PSIM phase residency, so it does not
extend `PSIM` lifetime for unrelated diagnostics.

Run directories:

```
Spark: tests/profile/si64/runs/setup-psim-isolated-spark-nstep1-20260601-140034
Spark: tests/profile/si64/runs/setup-psim-isolated-spark-nstep3-20260601-140123
Terok: tests/profile/si64/runs/setup-psim-isolated-terok-nstep1-20260601-140034
Terok: tests/profile/si64/runs/setup-psim-isolated-terok-nstep3-20260601-140123
```

| System | Case | NSTEPS | Wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | --- |
| Spark GB10 | `gpu_resident_stack` | 1 | 12.37 s | 1.3892 GB | yes |
| Spark GB10 | `gpu_resident_stack_setup_psim_host` | 1 | 12.44 s | 1.5102 GB | yes |
| Terok A40 | `gpu_resident_stack` | 1 | 18.06 s | 1.3892 GB | yes |
| Terok A40 | `gpu_resident_stack_setup_psim_host` | 1 | 15.47 s | 1.5102 GB | yes |
| Spark GB10 | `gpu_resident_stack` | 3 | 27.96 s | 2.8783 GB | n/a |
| Spark GB10 | `gpu_resident_stack_setup_psim_host` | 3 | 28.23 s | 2.9994 GB | n/a |
| Terok A40 | `gpu_resident_stack` | 3 | 29.54 s | 2.8783 GB | n/a |
| Terok A40 | `gpu_resident_stack_setup_psim_host` | 3 | 39.17 s | 2.9994 GB | n/a |

The copy-boundary evidence is the important part of this diagnostic. With setup
`PSIM` residency enabled, both Spark and Terok record
`ACC_COPY_SETUP_PSIM_IN` once and `ACC_PRESENT_PROP_PSIM` for all three
propagation calls in the NSTEPS=3 run. With the ablation enabled, the profile
falls back to `ACC_COPY_GRAM_PSIM_PSI_IN`, `ACC_COPY_GRAM_PSIM_PSI_OUT`, and one
`ACC_COPY_PROP_PSIM_IN` row. This removes one 0.1210 GB propagation input copy
from the stack default. The single-step wall time is noisy, especially on Terok,
but the three-step run confirms that the reduced copy boundary is not just an
accounting artifact.

## Previous Full Matrix Comparison

The table below is retained as the May 30 full-matrix comparison. The current
June 1 refresh matrix is listed in the next section.

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

## 2026-06-01 Refresh Matrix

This Spark C86C refresh used current `cp-paw-nvhpc` commit `85882ef` with
`TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`, one GPU rank for GPU cases,
and one-rank plus eight-rank CPU references. All listed cases completed and
matched the Si64 reference energy within tolerance.

| Suite | Case | Ranks | Wall time | Copy estimate | Interpretation |
| --- | --- | ---: | ---: | ---: | --- |
| GPU | `gpu_resident_hpsi_opsi` | 1 | 12.83 s | 1.96 GB | Best current case; broad PSI/HPSI/OPSI residency wins. |
| GPU | `gpu_resident` | 1 | 12.88 s | 2.69 GB | Same performance class with slightly more copy traffic. |
| GPU | `gpu_resident_hpsi` | 1 | 12.90 s | 2.20 GB | HPSI residency is correct and essentially tied. |
| GPU | `gpu_no_cufft` | 1 | 45.01 s | 96.55 GB | cuFFT is not the useful lever; cuBLAS still matters. |
| GPU | `cublas` | 1 | 49.71 s | 96.45 GB | Explicit cuBLAS alone is good, but residency is much better. |
| GPU | `gpu_no_cusolver` | 1 | 50.56 s | 104.99 GB | cuSOLVER is not decisive for this Si64 matrix. |
| GPU | `cusolver` | 1 | 70.26 s | 0.10 GB | Correct, but near CPU wall time without residency. |
| GPU | `gpu_off` | 1 | 70.27 s | 0.00 GB | Same GPU-capable binary with accelerators disabled. |
| GPU | `cufft` | 1 | 70.79 s | 0.00 GB | Native cuFFT alone does not help this workload. |
| GPU | `nvlamath` | 1 | 72.81 s | 0.00 GB | Builds and runs after the Cholesky symbol fix, but no speedup. |
| GPU | `gpu_no_cublas` | 1 | 73.31 s | 8.64 GB | Disabling cuBLAS removes most of the GPU benefit. |
| GPU | `nvblas` | 1 | 91.50 s | 0.00 GB | Interposition path remains too slow as a default. |
| GPU | `gpu_all` | 1 | 155.46 s | 96.55 GB | Combining all libraries is counterproductive here. |
| GPU | `cufftw` | 1 | 177.06 s | 0.00 GB | FFTW-compatible cuFFT wrapper is not attractive. |
| GPU | `gpu_all_off` | 1 | 178.07 s | 0.00 GB | All-library binary with accelerators disabled is diagnostic only. |
| CPU | `nvhpc_cpu` | 1 | 69.50 s | 0.00 GB | Best one-rank CPU reference. |
| CPU | `cpu` | 1 | 73.32 s | 0.00 GB | Plain one-rank CPU reference. |
| CPU | `nvhpc_cpu` | 8 | 166.88 s | 0.00 GB | Eight-rank CPU/NVHPC reference; slower wall time for this case. |
| CPU | `cpu` | 8 | 167.49 s | 0.00 GB | Eight-rank plain CPU reference; MPI overhead dominates. |

The decisive comparison is therefore `gpu_resident_hpsi_opsi` at 12.83 s
against `nvhpc_cpu` at 69.50 s on one rank and 166.88 s on eight ranks. This
strongly favors continued wavefunction/projector residency work over adding
more optional libraries to the default path.

## 2026-06-01 4fbe2cd Standard Refresh

After adding the direct combined HPSI/OPSI benchmark cases, Spark C86C was
rerun at current `cp-paw-nvhpc` commit `4fbe2cd` with `TEST=si64_bands`,
`EMPTY_BANDS=1024`, `NSTEPS=1`, one GPU rank for GPU cases, and one-rank plus
eight-rank CPU references. The default standard GPU list now includes the safe
HPSI/OPSI residency path and the full focused residency stack:

```
runs/si64_bands-nvhpc-standard-20260601-4fbe2cd-1024-nstep1
```

| Suite | Case | Ranks | Wall time | Copy estimate | Energy check | Interpretation |
| --- | --- | ---: | ---: | ---: | --- | --- |
| GPU | `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 1 | 12.21 s | 1.9940 GB | yes | Best 1024/1 case; HPSI, OPSI, PROJ, DENMAT energy, and off-site device-pack accumulation combine cleanly. |
| GPU | `gpu_resident_hpsi_opsi` | 1 | 12.77 s | 1.9558 GB | yes | Safest broad wavefunction-residency comparison; very close to the full stack. |
| GPU | `gpu_resident` | 1 | 12.87 s | 2.6872 GB | yes | Baseline resident path remains strong but moves more data. |
| GPU | `gpu_resident_hpsi` | 1 | 12.94 s | 2.1978 GB | yes | HPSI-only residency is correct but not the best 1024/1 point. |
| GPU | `gpu_resident_no_cusolver` | 1 | 15.66 s | 2.5916 GB | yes | cuSOLVER still helps the current resident matrix, but less than residency itself. |
| GPU | `gpu_resident_invbatch_off` | 1 | 16.22 s | 73.2764 GB | no | Negative-control case; disabling inversion batching changes the Si64 energy and must not be promoted. |
| GPU | `gpu_off` | 1 | 77.94 s | 0.0000 GB | yes | Same GPU-capable binary with accelerators disabled; near the plain CPU reference. |
| CPU | `nvhpc_cpu` | 1 | 75.23 s | 0.0000 GB | yes | Best one-rank CPU reference in this refresh. |
| CPU | `cpu` | 1 | 77.79 s | 0.0000 GB | yes | Plain one-rank CPU reference. |
| CPU | `nvhpc_cpu` | 8 | 392.07 s | 0.0000 GB | yes | NVHPC/NVPL improves the eight-rank CPU run, but MPI/setup overhead dominates this small smoke. |
| CPU | `cpu` | 8 | 672.52 s | 0.0000 GB | yes | Plain eight-rank CPU reference; not a useful resource-equivalent comparison for this case. |

The key comparison is now `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum`
at 12.21 s against `gpu_resident_hpsi_opsi` at 12.77 s and the one-rank NVHPC
CPU reference at 75.23 s. The focused full stack is energy-valid and finally
has a positive 1024/1 signal, but it should still be treated as a candidate for
the next larger case rather than an immediate production default. The 8-rank
CPU numbers are included for completeness; in this small smoke they mostly
measure MPI/setup overhead and should not be used to reject GPU offload.

## 2026-06-01 2048-Band Focus Validation

The follow-up focused run used current `cp-paw-nvhpc` commit `0a43a65` with
`TEST=si64_bands`, `EMPTY_BANDS=2048`, `NSTEPS=1`, and one GPU rank:

```
runs/si64_bands-focus-20260601-0a43a65-2048-nstep1-1r
```

| Case | Wall time | Copy estimate | Energy check | Interpretation |
| --- | ---: | ---: | --- | --- |
| `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 35.28 s | 4.0908 GB | yes | Best focused case; the 1024 gain is not just a small-case artifact. |
| `gpu_resident_hpsi_opsi` | 37.17 s | 3.9844 GB | yes | Safe broad residency comparison remains close but slower. |
| `gpu_resident` | 37.33 s | 5.3749 GB | yes | Baseline residency path is energy-valid but moves more data. |

This is enough evidence to add the shorter harness case `gpu_resident_stack`
and the runtime meta-keyword `CPPAW_GPU_RESIDENCY_STACK=1`. The keyword enables
the same focused HPSI/OPSI/PROJ/DENMAT/off-site device-pack accumulation stack
while leaving the individual `CPPAW_GPU_*` switches available for ablations.

The keyword was smoke-tested after adding the lowered stack thresholds:

| Run directory | Ranks | `gpu_resident_stack` | Long explicit stack | Energy check | Notes |
| --- | ---: | ---: | ---: | --- | --- |
| `residency-stack-keyword-20260601-fixed-512-nstep1-1r` | 1 | 6.22 s | 6.74 s | yes | Same energy and copy estimate as the long explicit case. |
| `residency-stack-keyword-20260601-fixed-512-nstep1-4r` | 4 | 16.77 s | 9.70 s | yes | Both cases use the same profile path, including `CUBLAS_ZGEMM_PROJ_RES`; this short four-rank smoke is run-order/noise sensitive. |
| `residency-stack-keyword-20260601-fixed-512-nstep1-4r-rev` | 4 | 33.65 s | 32.28 s | yes | Reversing the case order makes the wall times converge, confirming that the keyword is not missing the projection/off-site stack. |

## VPSI Boundary Harness

The next focused harness is `run_vpsi_boundary.sh`. It compares
`gpu_resident_hpsi`, `gpu_resident_hpsi_opsi`, and `gpu_resident_stack` around
the `WAVES_VPSI` producer boundary, then writes both the normal benchmark
summary and selected profile rows. `profile_copy_rows.py` now accepts
`--op-prefix`, `--op-regex`, `--include-zero`, and `--sort-by`, so the harness
can report copy/present/update rows separately from seconds-sorted
`PAW_VPSI_*` timing rows.

Spark C86C validation:

```
runs/vpsi-boundary-20260601-d3cd6cc-512-nstep1
runs/vpsi-boundary-20260601-7ad2625-smoke
```

| Suite | Case | Ranks | Wall time | Copy estimate | Energy check |
| --- | --- | ---: | ---: | ---: | --- |
| 512/NSTEPS=1 | `gpu_resident_hpsi` | 1 | 6.00 s | 1.2332 GB | yes |
| 512/NSTEPS=1 | `gpu_resident_hpsi_opsi` | 1 | 6.85 s | 1.0988 GB | yes |
| 512/NSTEPS=1 | `gpu_resident_stack` | 1 | 6.56 s | 1.1155 GB | yes |
| 512/NSTEPS=1 | `gpu_resident_hpsi` | 4 | 31.37 s | 1.5941 GB | yes |
| 512/NSTEPS=1 | `gpu_resident_hpsi_opsi` | 4 | 31.75 s | 1.4596 GB | yes |
| 512/NSTEPS=1 | `gpu_resident_stack` | 4 | 28.46 s | 1.5088 GB | yes |

The final one-case reporting smoke for `gpu_resident_stack` produced these
VPSI timing rows:

| Row | Seconds |
| --- | ---: |
| `PAW_VPSI_TOTAL` | 0.7235 |
| `PAW_VPSI_FFT_GTOR` | 0.3565 |
| `PAW_VPSI_FFT_RTOG` | 0.3565 |
| `PAW_VPSI_POT` | 0.0083 |
| `PAW_VPSI_KIN` | 0.0009 |

This confirms the design direction: optimizing the scalar real-space potential
or kinetic loops inside `WAVES_VPSI` will not move the needle for Si64. The next
actual implementation target should be GPU-resident FFT/RTOG or a broader
producer-side wavefunction region that removes the required `HPSI` refresh
after host-side FFT output.

## PSI0-To-PRINFO Residency

The next positive cross-boundary residency step keeps `PSI0` present after the
OPSI build until `PRINFO/WRITEPDOS` has consumed it for the PDOS projection.
`WAVES$SWITCH` deletes any still-resident old `PSI0` before swapping the
`PSI0`/`PSIM` pointers, so the OpenACC present table cannot carry stale data
into the next step. The opt-out diagnostic case is
`gpu_resident_stack_psi0_prinfo_host`, equivalent to
`CPPAW_GPU_RESIDENCY_STACK=1 CPPAW_GPU_PSI0_PRINFO_RESIDENCY=0`.

Validation runs used `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=1`,
`REPEATS=3`, one GPU rank:

```
runs/psi0-prinfo-spark-20260601-074852
runs/psi0-prinfo-terok-20260601-074852
```

| Machine | Case | Median wall time | Transfer estimate | Key marker | Energy check |
| --- | --- | ---: | ---: | --- | --- |
| Spark C86C | `gpu_resident_stack` | 12.11 s | 1.5099 GB | `ACC_PRESENT_PROJ_WRITEPDOS_PSI` | yes |
| Spark C86C | `gpu_resident_stack_psi0_prinfo_host` | 12.24 s | 1.6309 GB | `ACC_COPY_PROJ_WRITEPDOS_PSI_IN` | yes |
| Terok A40 | `gpu_resident_stack` | 12.77 s | 1.5099 GB | `ACC_PRESENT_PROJ_WRITEPDOS_PSI` | yes |
| Terok A40 | `gpu_resident_stack_psi0_prinfo_host` | 14.03 s | 1.6309 GB | `ACC_COPY_PROJ_WRITEPDOS_PSI_IN` | yes |

The copy reduction is deterministic: the stack removes exactly one
`THIS%PSI0` projection input copy for this 1024-band Si64 step, 0.1210 GB. Wall
time remains too noisy for a speedup claim, but this is the first clean
post-orthogonalization consumer reuse and it validates the broader residency
design Peter suggested: keep wavefunction data resident across multiple PAW
phases, with explicit lifecycle cleanup at pointer-swap boundaries.

## HPSI-To-Propagate Residency Diagnostic

The next targeted diagnostic keeps `HPSI` resident after the ETOT
expectation/Hamiltonian overlaps when the GPU `PSIM` propagation path is active.
This is controlled by `CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1` and the harness case
`gpu_resident_stack_hpsi_prop_psim_phase`. It is deliberately not enabled by the
plain stack keyword yet, because it changes the propagation path and still needs
the surrounding `PSIM` lifecycle to become more complete.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries with `TEST=si64_bands`,
`EMPTY_BANDS=1024`, `NSTEPS=1`, `REPEATS=3`, one GPU rank:

```
runs/hpsi-prop-spark-20260601-100735
runs/hpsi-prop-terok-20260601-101036
```

| Machine | Case | Median wall time | Transfer estimate | Key marker | Energy check |
| --- | --- | ---: | ---: | --- | --- |
| Spark C86C | `gpu_resident_stack` | 12.25 s | 1.5099 GB | baseline stack | yes |
| Spark C86C | `gpu_resident_stack_psim_phase` | 12.23 s | 1.6312 GB | `ACC_COPY_PROP_HPSI_IN` | yes |
| Spark C86C | `gpu_resident_stack_hpsi_prop_psim_phase` | 12.13 s | 1.5102 GB | `ACC_PRESENT_PROP_HPSI` | yes |
| Terok A40 | `gpu_resident_stack` | 13.75 s | 1.5099 GB | baseline stack | yes |
| Terok A40 | `gpu_resident_stack_psim_phase` | 12.33 s | 1.6312 GB | `ACC_COPY_PROP_HPSI_IN` | yes |
| Terok A40 | `gpu_resident_stack_hpsi_prop_psim_phase` | 12.59 s | 1.5102 GB | `ACC_PRESENT_PROP_HPSI` | yes |

The profile rows are identical on Spark and Terok: over three repeats,
`gpu_resident_stack_psim_phase` records `ACC_COPY_PROP_HPSI_IN` at 0.3631 GB,
while the new case records `ACC_PRESENT_PROP_HPSI` and no `PROP_HPSI` copy. This
turns the PSIM-phase path from a net-copy regression into a near-transfer-neutral
diagnostic. Wall time remains noisy, so the result is a validated residency
building block rather than a default promotion.

## PSIM-To-Switch Residency Diagnostic

The next cross-step diagnostic carries the resident `PSIM` allocation through
`WAVES$SWITCH` and marks it as the next step's resident `PSI0` after the host
pointers are swapped. This is controlled by
`CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1` and exposed in the harness as
`gpu_resident_stack_hpsi_prop_psim_switch`. It depends on the preceding
PSIM/HPSI propagation switches and remains opt-in because multi-step
wavefunction residency still needs more lifecycle coverage.

This change also replaces the old OpenACC directive delete in `WAVES$SWITCH`
with the OpenACC runtime `acc_delete` call for still-resident old `PSI0`/`PSIM`
arrays. Spark exposed a partially-present failure for the directive form at
`NSTEPS=2`; the runtime delete path completed on both Spark and Terok.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries with `TEST=si64_bands`,
`EMPTY_BANDS=1024`, `NSTEPS=2`, `EXPECTED_ENERGY=269.022536`, `REPEATS=3`, and
one GPU rank:

```
runs/psim-switch-accdelete-spark-20260601-104028
runs/psim-switch-accdelete-terok-20260601-104027
```

| Machine | Case | Median wall time | Transfer estimate | Key marker | Energy check |
| --- | --- | ---: | ---: | --- | --- |
| Spark C86C | `gpu_resident_stack_hpsi_prop_psim_phase` | 20.84 s | 2.5352 GB | baseline two-step PSIM phase | yes |
| Spark C86C | `gpu_resident_stack_hpsi_prop_psim_switch` | 20.61 s | 2.4142 GB | `ACC_PRESENT_SWITCH_PSI0_PSIM` | yes |
| Terok A40 | `gpu_resident_stack_hpsi_prop_psim_phase` | 20.96 s | 2.5352 GB | baseline two-step PSIM phase | yes |
| Terok A40 | `gpu_resident_stack_hpsi_prop_psim_switch` | 20.25 s | 2.4142 GB | `ACC_PRESENT_SWITCH_PSI0_PSIM` | yes |

The deterministic effect is modest but clean: the switch case removes the
force-side `ACC_COPY_FORCE_PSI0_IN` row, 0.3631 GB over three repeats, or
0.1210 GB per two-step run. Wall time moves in the right direction on both
machines, but the important result is correctness plus an explicit cross-step
wavefunction lifecycle hook for the next residency pass.

After shortening the profile marker to avoid CSV truncation, one-repeat marker
smokes on both machines produced `ACC_PRESENT_SWITCH_PSI0_PSIM` with `ok=yes`
and the same 2.4142 GB transfer estimate.

## Resident Wavefunction Lifecycle Dimensions

The follow-up refactor stores explicit resident dimensions for `PSI0`, `PSIM`,
and `HPSI` in each `WVSET_TYPE` and routes the common OpenACC lifecycle changes
through mark/clear/delete helpers. This keeps `WAVES$SWITCH` and later cleanup
sites from recomputing resident array sizes from the current `GSET`/`PROJ`
state, which was exactly the fragile boundary exposed by the earlier
partially-present Spark failure.

This is intended as a robustness/enabling step rather than a new performance
claim. It preserves the existing transfer accounting while making the next
cross-step residency extensions less dependent on pointer-swap timing.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries on Spark C86C and Terok.
The serial runs used `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=2`,
`EXPECTED_ENERGY=269.022536`, `REPEATS=3`, and one GPU rank:

```
runs/accdims-switch-spark-20260601-105947
runs/accdims-switch-terok-20260601-105946
```

| Machine | Case | Median wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | --- |
| Spark C86C | `gpu_resident_stack_hpsi_prop_psim_phase` | 20.82 s | 2.5352 GB | yes |
| Spark C86C | `gpu_resident_stack_hpsi_prop_psim_switch` | 20.84 s | 2.4142 GB | yes |
| Terok A40 | `gpu_resident_stack_hpsi_prop_psim_phase` | 24.54 s | 2.5352 GB | yes |
| Terok A40 | `gpu_resident_stack_hpsi_prop_psim_switch` | 23.50 s | 2.4142 GB | yes |

Parallel smoke checks used `EMPTY_BANDS=512`, `NSTEPS=2`, `RANKS=4`, and the
switch-residency case:

| Machine | Run directory | Wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | --- |
| Spark C86C | `runs/accdims-switch-parallel-smoke-spark-20260601-110234` | 15.04 s | 2.0691 GB | yes |
| Terok A40 | `runs/accdims-switch-parallel-smoke-terok-20260601-110233` | 41.20 s | 2.0691 GB | yes |

The transfer estimates remain identical to the pre-refactor switch run:
2.5352 GB for the phase-only diagnostic and 2.4142 GB for the switch-resident
case. That is the desired result: this patch tightens lifecycle bookkeeping
without changing the numerical path or promoting a new default.

## Current VPSI/cuFFT Refresh

After adding explicit wavefunction residency dimensions, the VPSI boundary
harness was rerun to check whether the latest stack changes alter the cuFFT
decision. The focused comparison used `gpu_resident_stack`,
`gpu_resident_stack_cufft`, and `gpu_resident_stack_cufft_force` with
`NSTEPS=2`, `REPEATS=1`, one GPU rank, and energy checking on Spark C86C and
Terok:

```
runs/vpsi-cufft-refresh-spark-20260601-110633
runs/vpsi-cufft-refresh-terok-20260601-110633
```

| Machine | Empty bands | Case | Wall time | Transfer estimate | `PAW_VPSI_TOTAL` | Energy check |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| Spark C86C | 512 | `gpu_resident_stack` | 10.79 s | 1.4366 GB | 1.4492 s | yes |
| Spark C86C | 512 | `gpu_resident_stack_cufft` | 10.14 s | 1.4366 GB | 1.4750 s | yes |
| Spark C86C | 512 | `gpu_resident_stack_cufft_force` | 13.37 s | 10.9731 GB | 3.6147 s | yes |
| Spark C86C | 1024 | `gpu_resident_stack` | 20.92 s | 2.5346 GB | 2.6477 s | yes |
| Spark C86C | 1024 | `gpu_resident_stack_cufft` | 20.77 s | 2.5346 GB | 2.6498 s | yes |
| Spark C86C | 1024 | `gpu_resident_stack_cufft_force` | 26.31 s | 19.6164 GB | 6.3722 s | yes |
| Terok A40 | 512 | `gpu_resident_stack` | 10.80 s | 1.4366 GB | 1.5866 s | yes |
| Terok A40 | 512 | `gpu_resident_stack_cufft` | 10.33 s | 1.4366 GB | 1.4988 s | yes |
| Terok A40 | 512 | `gpu_resident_stack_cufft_force` | 19.62 s | 10.9731 GB | 7.8295 s | yes |
| Terok A40 | 1024 | `gpu_resident_stack` | 27.51 s | 2.5346 GB | 2.6961 s | yes |
| Terok A40 | 1024 | `gpu_resident_stack_cufft` | 20.79 s | 2.5346 GB | 2.5499 s | yes |
| Terok A40 | 1024 | `gpu_resident_stack_cufft_force` | 38.83 s | 19.6164 GB | 13.0076 s | yes |

Interpretation: threshold-gated native cuFFT remains safe as a diagnostic and
can be neutral to slightly favorable in this one-repeat refresh, but the
force-all mode is decisively worse. It increases the estimated transfer volume
by roughly 9.5 GB at 512 empty bands and 17.1 GB at 1024 empty bands, and it
more than doubles `PAW_VPSI_TOTAL` on both machines. The next useful FFT work is
therefore not broader forced cuFFT activation inside the current host-oriented
3D FFT path. It should be either a genuinely device-resident FFT/RTOG pipeline
or a narrower effort to keep the wavefunction buffers resident around the
current FFT producer/consumer boundary.

## Lazy Projector Scratch Allocation

The resident `PRO` cache path no longer allocates and fills host-side
`GVEC`/`PRO`/`EIGR` scratch arrays before entering the cached GPU projection and
addproduct paths. Those arrays are still allocated lazily for the host-expanded
fallback and for `CPPAW_GPU_PRO_EXPANSION=0`, so the ablation path remains
available.

Validation rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` on Spark C86C and Terok, then checked
both the recommended cache path and the host-PRO ablation:

```
runs/lazy-scratch-spark-20260601-1118
runs/lazy-scratch-terok-20260601-1118
runs/lazy-scratch-1024-spark-20260601-1121
runs/lazy-scratch-1024-terok-20260601-1121
runs/lazy-scratch-parallel-spark-20260601-1118
runs/lazy-scratch-parallel-terok-20260601-1120
```

| Machine | Empty bands | Ranks | Case | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| Spark C86C | 512 | 1 | `gpu_resident_stack` | 6.63 s | 0.8465 GB | yes |
| Spark C86C | 512 | 1 | `gpu_resident_pro_host` | 6.80 s | 2.9023 GB | yes |
| Terok A40 | 512 | 1 | `gpu_resident_stack` | 6.72 s | 0.8465 GB | yes |
| Terok A40 | 512 | 1 | `gpu_resident_pro_host` | 7.59 s | 2.9023 GB | yes |
| Spark C86C | 1024 | 1 | `gpu_resident_stack` | 12.12 s | 1.5099 GB | yes |
| Spark C86C | 1024 | 1 | `gpu_resident_pro_host` | 13.01 s | 4.0858 GB | yes |
| Terok A40 | 1024 | 1 | `gpu_resident_stack` | 11.99 s | 1.5099 GB | yes |
| Terok A40 | 1024 | 1 | `gpu_resident_pro_host` | 14.36 s | 4.0858 GB | yes |
| Spark C86C | 512 | 4 | `gpu_resident_stack` | 9.63 s | 1.2398 GB | yes |
| Terok A40 | 512 | 4 | `gpu_resident_stack` | 8.14 s | 1.2398 GB | yes |

Interpretation: this is a cleanup and robustness step for the recommended
resident-cache stack, not a new library or a claimed standalone speedup. It
removes unused host work from the hot projection/addproduct setup while keeping
the `gpu_resident_pro_host` fallback exercised and correct.

## PSIM Stack-Default Probe

The existing opt-in PSIM/HPSI propagation and switch residency was tested as a
candidate for promotion into the broad `CPPAW_GPU_RESIDENCY_STACK` default. A
temporary build enabled:

```
CPPAW_GPU_PSIM_PROPAGATE=1
CPPAW_GPU_PSIM_PHASE_RESIDENCY=1
CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1
CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1
```

inside the stack keyword, while `gpu_resident_stack_legacy_psim` disabled those
switches again as the old-stack ablation.

Validation/probe runs:

```
runs/psim-promote-probe-spark-20260601-1130
runs/psim-promote-probe-terok-20260601-1130
runs/psim-stack-default-spark-20260601-1140
runs/psim-stack-default-terok-20260601-1140
runs/psim-stack-default-parallel-spark-20260601-1145
runs/psim-stack-default-parallel-terok-20260601-1145
runs/psim-stack-default-parallel-terok-repeat-20260601-1150
```

| Machine | Empty bands | NSTEPS | Ranks | Case | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| Spark C86C | 1024 | 2 | 1 | promoted stack | 21.29 s | 2.4142 GB | yes |
| Spark C86C | 1024 | 2 | 1 | legacy stack | 21.18 s | 2.5346 GB | yes |
| Terok A40 | 1024 | 2 | 1 | promoted stack | 20.77 s | 2.4142 GB | yes |
| Terok A40 | 1024 | 2 | 1 | legacy stack | 29.17 s | 2.5346 GB | yes |
| Spark C86C | 512 | 2 | 4 | promoted stack | 15.14 s | 2.0691 GB | yes |
| Spark C86C | 512 | 2 | 4 | legacy stack | 15.21 s | 2.1357 GB | yes |
| Terok A40 | 512 | 2 | 4 | promoted stack | 41.44 s / 37.94 s | 2.0691 GB | yes |
| Terok A40 | 512 | 2 | 4 | legacy stack | 12.33 s / 23.15 s | 2.1357 GB | yes |

Interpretation at the time: the opt-in path was useful for the intended
resource split of one MPI rank per GPU, but this pre-cleanup promotion probe
showed a Terok 4-rank shared-GPU slowdown. The later resident-dimension and
lifecycle cleanup removed that instability; see the next section for the
current default decision.

## PSIM Switch In Stack Default

After the resident wavefunction dimension tracking and cleanup helpers were in
place, the PSIM propagation/switch stack was promoted again and retested as the
plain `gpu_resident_stack` default. The default now enables:

```
CPPAW_GPU_PSIM_PROPAGATE=1
CPPAW_GPU_PSIM_PHASE_RESIDENCY=1
CPPAW_GPU_HPSI_PROPAGATE_RESIDENCY=1
CPPAW_GPU_PSIM_SWITCH_RESIDENCY=1
```

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries on Spark C86C and Terok:

```
runs/stack-default-psim-spark-512-nstep2-1r-20260601-120645
runs/stack-default-psim-terok-512-nstep2-1r-20260601-120645
runs/stack-default-psim-spark-512-nstep2-4r-20260601-120826
runs/stack-default-psim-terok-512-nstep2-4r-20260601-120826
runs/stack-default-psim-spark-512-nstep1-1r-20260601-120855
runs/stack-default-psim-terok-512-nstep1-1r-20260601-120855
```

| Machine | Empty bands | NSTEPS | Ranks | Case | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| Spark C86C | 512 | 2 | 1 | `gpu_resident_stack` avg | 9.73 s | 1.3700 GB | yes |
| Spark C86C | 512 | 2 | 1 | explicit switch avg | 9.62 s | 1.3700 GB | yes |
| Terok A40 | 512 | 2 | 1 | `gpu_resident_stack` avg | 11.07 s | 1.3700 GB | yes |
| Terok A40 | 512 | 2 | 1 | explicit switch avg | 11.06 s | 1.3700 GB | yes |
| Spark C86C | 512 | 2 | 4 | `gpu_resident_stack` | 15.07 s | 2.0691 GB | yes |
| Terok A40 | 512 | 2 | 4 | `gpu_resident_stack` | 12.03 s | 2.0691 GB | yes |
| Spark C86C | 512 | 1 | 1 | `gpu_resident_stack` | 6.65 s | 0.8468 GB | yes |
| Terok A40 | 512 | 1 | 1 | `gpu_resident_stack` | 6.61 s | 0.8468 GB | yes |

Profile checks confirm that `ACC_COPY_ORTHO_PSIM_IN` disappears from the plain
`gpu_resident_stack` profile, while `ACC_COPY_ORTHO_PSIM_OUT` remains because
the host still needs the updated wavefunction after orthogonalization. Compared
with the old stack at 512 empty bands and `NSTEPS=2`, this removes 0.0666 GB of
deterministic transfer on both machines. The explicit
`gpu_resident_stack_hpsi_prop_psim_switch` case is now equivalent to the stack
default and remains useful only as a compatibility/debug spelling.

## FFT/VPSI Collector Fields

The benchmark tooling now exposes the FFT/VPSI structure directly in the main
`benchmark.tsv` output instead of requiring a manual profile-row query. New
columns are:

```
pw_fft_gtor_s
pw_fft_rtog_s
vpsi_s
vpsi_gtor_s
vpsi_rtog_s
```

The focused VPSI harness also includes `PW_FFT_*` rows in its FFT-phase tables.
Re-summarizing the existing 1024/1 lazy-scratch profiles gives:

| Machine | Case | `fft_s` | `pw_fft_gtor_s` | `pw_fft_rtog_s` | `vpsi_s` | `vpsi_gtor_s` | `vpsi_rtog_s` |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Spark C86C | `gpu_resident_stack` | 3.6887 s | 1.4310 s | 0.6335 s | 1.2877 s | 0.6415 s | 0.6279 s |
| Spark C86C | `gpu_resident_pro_host` | 3.8184 s | 1.4800 s | 0.6640 s | 1.3339 s | 0.6517 s | 0.6583 s |
| Terok A40 | `gpu_resident_stack` | 3.2424 s | 1.4389 s | 0.6048 s | 1.2514 s | 0.6135 s | 0.5964 s |
| Terok A40 | `gpu_resident_pro_host` | 3.6762 s | 1.6118 s | 0.8601 s | 1.8034 s | 0.8833 s | 0.8513 s |

This does not change physics or runtime behavior. It makes future night runs
more diagnostic: a useful FFT/RTOG change should reduce `vpsi_*` and the
corresponding `PW_FFT_*` columns, not only move time between broad `PAW_*`
envelopes.

## Current Conclusions

1. Use the residency profile path as the recommended NVHPC GPU profiling path:
   `nvhpc_gpu_acc_residency_profile` and
   `nvhpc_gpu_acc_residency_profile_parallel`.

2. The main win is device residency around wavefunction-heavy regions plus
   explicit cuBLAS. In the latest 1024/1 standard run, the focused full-stack
   case is 12.21 s and `gpu_resident_hpsi_opsi` is 12.77 s versus 75.23 s for
   the one-rank NVHPC CPU reference. The 2048/1 focus run keeps the same order:
   35.28 s for the focused stack versus 37.17 s for HPSI/OPSI residency.

3. Do not make all optional NVIDIA libraries active by default. The
   `gpu_all*`, `cufftw`, `nvblas`, and `nvlamath` cases are valuable diagnostics
   but are slower for this workload. The focused combined residency stack is
   different from broad all-library activation and is now exposed as
   `CPPAW_GPU_RESIDENCY_STACK=1` for focused benchmark runs.

4. Keep native cuFFT and cuSOLVER threshold-gated. The Si64 result does not
   justify aggressive defaults for either one: native `cufft` is 70.79 s and
   `cusolver` is 70.26 s in the refresh matrix. Larger generalized eigensolver
   cases may still change the cuSOLVER decision.

5. The first projector follow-up is now in place: resident `PRO` is built once,
   cached on the GPU and reused by projection plus eligible addproduct calls.
   Keep that full path as part of the residency default because it reduces copy
   estimates substantially, while the latest benchmark shows that the later
   PSI/HPSI/OPSI residency work is the dominant Spark speedup. The latest
   lazy-scratch cleanup keeps this path leaner by avoiding unused host
   projector scratch allocation in the cached resident `PRO` route.

6. The one-center overlap GPU-pack path removes the previous large
   `WAVES_1COVERLAP` bottleneck. The new `PSI0`-to-`PRINFO` and
   `PSIM`-to-switch residency diagnostics each remove another deterministic
   0.1210 GB copy in their focused runs, but wall time is still neutral/noisy.
   The PSIM/HPSI propagation path remains opt-in rather than part of
   `CPPAW_GPU_RESIDENCY_STACK=1` because Terok 4-rank shared-GPU validation
   shows a repeatable wall-time regression when it is promoted. The next useful
   default-candidate work should target the remaining host-side FFT/RTOG and
   producer-side HPSI/projector boundaries, not broader default activation of
   cuFFT/cuFFTW/NVLAMATH/NVBLAS.

7. Keep forced cuFFT out of the recommended stack. The current VPSI/cuFFT
   refresh confirms that threshold-gated cuFFT is harmless, but force-all cuFFT
   is slower and transfer-heavy on both Spark and Terok. A useful FFT follow-up
   needs device-resident dataflow, not just more `LIB$FFTC8` calls routed
   through cuFFT. The benchmark tooling now reports `vpsi_*` and `PW_FFT_*`
   timing columns directly so those changes can be evaluated from the main TSV.

8. Keep `gpu_resident_invbatch_off` as a negative-control diagnostic only. In
   the 4fbe2cd refresh it produced 302.773536 Ha instead of 302.280854 Ha and
   therefore failed the energy guard.

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
8-rank CPU reference. The default case list now includes the current HPSI,
OPSI, and stacked residency candidates, plus the older addproduct/projector
host-side ablations for continuity. Add `RUN_CPU_REFERENCES=yes` when CPU
reference numbers are explicitly needed. The wrapper delegates to
`run_nvhpc_standard.sh`, so it writes the same combined benchmark, comparison,
transfer-row, and present-row reports.

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

## ADDPRO Fallback Copy Accounting

The follow-up records `ACC_COPY_ADDPRO_<ctx>_PSI_IO` before the
host-expansion/addproduct fallback enters its OpenACC `COPY(PSI)` region. This
puts the fallback path on the same semantic accounting basis as the resident
projector-cache path. Per-atom `PRO`/`PROPSI` transfers in the fallback path
remain in the generic cuBLAS `ZGEMM_NN` copy rows to avoid double counting.

Spark C86C validation:

```
runs/addpro-fallback-copyacct-20260531-512-1r
runs/addpro-fallback-copyacct-20260531-512-4r
runs/addpro-fallback-copyacct-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.39 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_addpro_hpsi_host` | 512 | 1 | 6.46 s | 1.6786 GB | 0.000000401 Ha |
| `gpu_resident_addpro_opsi_host` | 512 | 1 | 7.25 s | 1.6786 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 1 | 7.06 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_opsi_addpro_host` | 512 | 1 | 7.02 s | 1.6786 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.50 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_addpro_hpsi_host` | 512 | 4 | 9.61 s | 2.0443 GB | 0.000000401 Ha |
| `gpu_resident_addpro_opsi_host` | 512 | 4 | 9.61 s | 2.0443 GB | 0.000000401 Ha |
| `gpu_resident_opsi` | 512 | 4 | 9.50 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_opsi_addpro_host` | 512 | 4 | 9.39 s | 2.0443 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 39.44 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_addpro_hpsi_host` | 2048 | 1 | 40.02 s | 5.5498 GB | 0.000000407 Ha |
| `gpu_resident_addpro_opsi_host` | 2048 | 1 | 40.99 s | 5.5498 GB | 0.000000407 Ha |
| `gpu_resident_opsi` | 2048 | 1 | 40.09 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_opsi_addpro_host` | 2048 | 1 | 38.97 s | 5.5498 GB | 0.000000407 Ha |

All runs are energy-valid. The previous 2048-band copy decrease for the
cache-disabled contexts is now gone; those contexts show a higher copy estimate,
as expected when the fallback `PSI` input/output region is counted explicitly.
The context switches remain diagnostics, while the resident ADDPRO cache stays
the recommended default path.

## HPSI Residency Diagnostic

The next opt-in diagnostic keeps the Hamiltonian-side `HPSI` wavefunction
resident after `WAVES_ADDPRO` and reuses it for the immediate expectation and
full-Hamiltonian overlaps:

- `CPPAW_GPU_HPSI_RESIDENCY=1`
- benchmark case: `gpu_resident_hpsi`

Spark C86C validation:

```
runs/hpsi-residency-20260531-512-1r
runs/hpsi-residency-20260531-2048-1r
runs/hpsi-residency-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.62 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 1 | 6.29 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_resident` | 2048 | 1 | 39.89 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 40.29 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident` | 512 | 4 | 9.40 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.61 s | 1.7284 GB | 0.000000401 Ha |

All runs are energy-valid. The diagnostic consistently lowers the semantic copy
estimate and records `ACC_PRESENT_EXPECT_HPSI` and `ACC_PRESENT_HAMILTON_HPSI`
for the reused wavefunction. Wall time improves for the small 1-rank smoke, but
is neutral to slightly slower for the larger and 4-rank checks, so the switch
remains opt-in rather than a new default.

## HPSI + OPSI Benchmark Keyword

The harness follow-up adds the explicit combined keyword
`gpu_resident_hpsi_opsi` and records inherited `CPPAW_*` accelerator switches in
the benchmark `env` column. This avoids ambiguous runs where an extra shell
environment switch affects the executable but is not visible in the TSV.

Spark C86C validation:

```
runs/hpsi-opsi-keyword-20260531-512-1r
runs/hpsi-opsi-keyword-20260531-2048-1r
runs/hpsi-opsi-keyword-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_opsi` | 512 | 1 | 7.21 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_opsi` | 2048 | 1 | 39.75 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 4 | 9.65 s | 1.7284 GB | 0.000000401 Ha |

The combined keyword is energy-valid and reproduces the HPSI copy reduction,
but OPSI residency does not add another copy reduction for these Si64 smoke
cases. Keep it as a reproducible diagnostic combination, not a new default.

## One-Center Overlap Copy Split

The one-center overlap cuBLAS path now replaces the aggregate
`ACC_COPY_CUBLAS_1COV` estimate with separate transfer rows for packed
projector input and overlap-matrix output. The change is accounting-only; the
kernel path and total copy estimate are unchanged within rounding.

Spark C86C validation:

```
runs/1cov-split-20260531-2048-1r
runs/1cov-split-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | `1COV_PROJ_IN` | `1COV_CMAT_OUT` | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 2048 | 1 | 40.48 s | 4.8988 GB | 0.2173 GB | 0.1894 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.78 s | 1.7284 GB | 0.0160 GB | 0.0164 GB | 0.000000401 Ha |

For the large one-rank smoke, 1COV still contributes about 0.4067 GB in total,
but it is not the leading copy source after HPSI residency. The larger remaining
rows are still the wavefunction input/output regions around Gram,
orthogonalization, and ADDPRO, so the next performance work should target
those repeated wavefunction edges before optimizing 1COV matrix-output copies.

## Wavefunction IO Copy Split

The generic 3D wavefunction IO profiler now records separate input and output
rows for arrays that are updated in an OpenACC data region. This replaces
aggregate `*_PSI_IO` rows in Gram-Schmidt, the orthogonalization outer region,
`WAVES_ADDPRO`, and `WAVES_ADDOPSI` without changing the actual data regions or
the total semantic copy estimate.

Spark C86C validation:

```
runs/wave-io-split-20260531-2048-1r
runs/wave-io-split-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Example split rows | Energy delta |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `gpu_resident_hpsi` | 2048 | 1 | 39.22 s | 4.8988 GB | `GRAM_PSI0_PSI_IN/OUT`, `GRAM_PSIM_PSI_IN/OUT`, `ORTHO_PSIM_IN/OUT`, `ADDPRO_OPSI_PSI_IN/OUT` each 0.2286 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.72 s | 1.7284 GB | same rows each 0.0168 GB | 0.000000401 Ha |

No old `ACC_COPY_*_PSI_IO` rows remain in the checked profiles. The split shows
that the next meaningful optimization is not a single output-only writeback:
for the large Si64 smoke, Gram `PSI0`, Gram `PSIM`, orthogonalization `PSIM`,
and ADDPRO `OPSI` each still have equal input and output sides. A future
residency experiment should therefore try to keep one of these wavefunctions
resident across the producing and consuming phases, rather than only suppressing
copy-out accounting.

## PSIM Propagation Diagnostic

The PSIM propagation follow-up adds the opt-in switch
`CPPAW_GPU_PSIM_PROPAGATE=1` and benchmark cases `gpu_psim_propagate` and
`gpu_hpsi_psim_propagate`. This path runs the final
`PSIM = A*PSI0 + B*PSIM + C*HPSI` update of `WAVES$PROPAGATE` on the GPU and
copies `PSIM` back before the next orthogonalization block. The older
environment aliases `CPPAW_GPU_PSIM_RESIDENCY` and
`CPPAW_CUBLAS_ACC_PSIM_RESIDENCY` remain accepted for compatibility.

Spark C86C validation:

```
runs/psim-propagate-subroutine-20260531-2048-1r
runs/psim-propagate-final-20260531-2048-1r
runs/psim-propagate-final-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 2048 | 1 | 39.20 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_psim_propagate` | 2048 | 1 | 39.39 s | 6.2897 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 39.50 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_hpsi_psim_propagate` | 2048 | 1 | 39.45 s | 5.8136 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 27.28 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 4 | 25.39 s | 1.9977 GB | 0.000000401 Ha |

All runs are energy-valid. A direct OpenACC kernel on the `THIS%PSIM` component
was energy-invalid for this case; moving the update into a standalone
dummy-array subroutine fixed correctness. The diagnostic is therefore useful as
a regression harness for future residency work, but it is not a new default:
the current implementation adds `PSI0`, `PSIM`, `HPSI`, coefficient input
copies and a `PSIM` copy-out. The performance-positive path is still broader
projector, PRO, and wavefunction residency that lets later phases consume the
GPU-resident data instead of copying it back immediately.

## PSIM Present-Or-Copy Accounting

The follow-up changes the PSIM propagation data region from unconditional
`copy/copyin` clauses to `present_or_copy` for `PSIM` and `present_or_copyin`
for `PSI0`/`HPSI`. The profile bookkeeping now uses the existing in/out helper
for `PSIM`, so a future broader resident region records `ACC_PRESENT_PROP_PSIM`
instead of an unconditional output copy row.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/psim-present-copy-20260601-512-1r
runs/psim-present-copy-20260601-512-4r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 6.36 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 1 | 7.35 s | 1.7730 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 1 | 7.27 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 1 | 7.41 s | 1.6369 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.18 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 4 | 9.31 s | 2.1387 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.65 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 4 | 9.66 s | 1.9977 GB | 0.000000401 Ha |

The current Si64 path still records `ACC_COPY_PROP_PSIM_IN` and
`ACC_COPY_PROP_PSIM_OUT`, which proves that no enclosing phase keeps `PSIM`
present into `WAVES$PROPAGATE` yet. The code path is nevertheless now ready for
that next step: once orthogonalization/projection leaves `PSIM` resident, the
propagation kernel can reuse it without another local data-region rewrite.

## PSIM Focus Harness

The follow-up adds `tests/profile/si64/run_psim_focus.sh` so the PSIM/HPSI
propagation comparison is reproducible instead of being an ad-hoc case list. The
default sweep compares `gpu_resident`, `gpu_psim_propagate`,
`gpu_resident_hpsi`, and `gpu_hpsi_psim_propagate` at `NSTEPS=1`,
`EMPTY_BANDS=512`, first with one GPU rank and then with four ranks sharing the
same GPU. Set `RUN_LARGE_GPU=yes` to add the 2048-band one-rank check.

Spark C86C validation:

```
runs/psim-focus-harness-20260601-512-1r
runs/psim-focus-harness-20260601-512-4r
runs/psim-focus-harness-20260601-combined.tsv
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.13 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 1 | 6.20 s | 1.7730 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 1 | 6.38 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 1 | 6.17 s | 1.6369 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.29 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 4 | 9.19 s | 2.1387 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.77 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 4 | 9.56 s | 1.9977 GB | 0.000000401 Ha |

The PSIM-propagation cases can look favorable in individual 512-band smokes, but
they still copy more data than the corresponding no-PSIM-offload paths. Treat
this harness as the guardrail for the broader future step: keeping `PSIM`
resident across the orthogonalization/propagation boundary instead of copying it
back immediately.

## PSIM Cross-Phase Residency

The cross-phase diagnostic adds `CPPAW_GPU_PSIM_PHASE_RESIDENCY=1` and harness
cases `gpu_resident_psim_phase` and `gpu_resident_hpsi_psim_phase`. The path is
still opt-in. It relies on the current timestep ordering, where
`WAVES$PROPAGATE()` is followed immediately by `WAVES$ORTHOGONALIZE()`, so the
updated `PSIM` can remain device-resident until the orthogonalization block
finishes and copies the orthogonalized wavefunction back.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/psim-phase-residency-20260601-512-1r
runs/psim-phase-residency-20260601-512-4r
runs/psim-phase-residency-20260601-combined.tsv
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 1 | 7.08 s | 1.5037 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 1 | 7.11 s | 1.7730 GB | 0.000000401 Ha |
| `gpu_resident_psim_phase` | 512 | 1 | 6.96 s | 1.6385 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 1 | 6.47 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 1 | 7.72 s | 1.6369 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_psim_phase` | 512 | 1 | 7.50 s | 1.5024 GB | 0.000000401 Ha |
| `gpu_resident` | 512 | 4 | 9.24 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_psim_propagate` | 512 | 4 | 9.34 s | 2.1387 GB | 0.000000401 Ha |
| `gpu_resident_psim_phase` | 512 | 4 | 9.11 s | 2.0042 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.68 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 4 | 9.89 s | 1.9977 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_psim_phase` | 512 | 4 | 9.57 s | 1.8632 GB | 0.000000401 Ha |

Profile rows show the intended data-lifetime change. The non-phase PSIM cases
record `ACC_COPY_PROP_PSIM_IN` and `ACC_COPY_PROP_PSIM_OUT`. The phase cases
record `ACC_COPY_PROP_PSIM_IN`, then `ACC_PRESENT_ORTHO_PSIM`, and finally
`ACC_COPY_ORTHO_PSIM_OUT`. This removes one immediate propagation copy-out and
lets orthogonalization consume the propagated device array directly. The copy
estimate drops by about 0.1345 GB for the 512-band one-rank case and about
0.1345 GB per four-rank run in the shared-GPU smoke. Timings are neutral to
slightly favorable in this small case, so the next meaningful test is a larger
band run before considering promotion beyond diagnostic status.

The 2048-band one-rank follow-up keeps the same conclusion but with more useful
kernel sizes:

```
runs/psim-phase-residency-20260601-2048-2048-1r
runs/psim-phase-residency-20260601-2048-combined.tsv
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 2048 | 1 | 37.61 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_psim_propagate` | 2048 | 1 | 38.26 s | 6.2897 GB | 0.000000407 Ha |
| `gpu_resident_psim_phase` | 2048 | 1 | 37.75 s | 5.8325 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 37.87 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_hpsi_psim_propagate` | 2048 | 1 | 38.26 s | 5.8136 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_psim_phase` | 2048 | 1 | 39.18 s | 5.3563 GB | 0.000000407 Ha |

Cross-phase residency removes about 0.4572 GB from the plain PSIM propagation
path at 2048 bands, and it makes the PSIM path nearly neutral against
`gpu_resident` in wall time. It still copies more than the no-PSIM-propagation
baseline because `PSI0`, `HPSI`, and the propagation coefficients are not yet
resident across the same boundary. The HPSI-plus-PSIM phase combination reduces
copy volume too, but this single run was slower, so it remains a diagnostic.

## DENMAT Profiling Split

The one-center density-matrix follow-up splits the previous coarse
`PAW_ETOT_DENMAT` envelope into `PAW_DENMAT_*` rows for occupation setup, site
setup, inner density/energy loops, accumulation, MPI combine, and spin
conversion. It also adds coarse `PAW_OFFDEN_*` rows for off-site density-matrix
bookkeeping. This is instrumentation only; it does not change the accelerator
data path.

Spark C86C validation:

```
runs/denmat-profiling-final-20260531-2048-1r
runs/denmat-profiling-final-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 2048 | 1 | 41.70 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.85 s | 1.7284 GB | 0.000000401 Ha |

| Profile row | 2048 bands, 1 rank | 512 bands, 4 ranks, rank 1 |
| --- | ---: | ---: |
| `PAW_ETOT_DENMAT` | 4.1846 s | 0.2378 s |
| `PAW_DENMAT_SITE_KERNEL` | 3.4608 s | 0.0827 s |
| `PAW_DENMAT_ENERGY_LOOP` | 3.0668 s | 0.0798 s |
| `PAW_DENMAT_DENSITY_LOOP` | 0.0131 s | 0.0023 s |
| `PAW_OFFDEN_SUM` | 0.7213 s | 0.1534 s |

All runs are energy-valid. The split shows that the large one-rank Si64 case is
not dominated by projector setup or accumulation inside `WAVES$DENMAT`; it is
dominated by the inner `WAVES_DENMAT` energy/Lambda loop. That makes a future
GPU or BLAS-style rewrite of the one-center energy-density contraction a better
candidate than further micro-optimizing DENMAT host setup. In the 4-rank smoke,
the local DENMAT kernel is much smaller and off-site summation is the larger
remaining subpiece, so a parallel follow-up should keep MPI/off-site behavior in
view.

## DENMAT Energy OpenACC Diagnostic

The follow-up adds an opt-in OpenACC path for the time-inversion
`WAVES_DENMAT` energy/Lambda contraction. The runtime switch is
`CPPAW_GPU_DENMAT_ENERGY=1` and the benchmark cases are
`gpu_resident_denmat_energy` and `gpu_resident_hpsi_denmat_energy`. The
implementation is two-stage: first build the per-band `FUNC` contraction on the
GPU, then form the small one-center output block. The path is disabled by
default because it still copies the full Lambda block per site.

Spark C86C validation:

```
runs/denmat-energy-acc-v2-20260531-512-1r
runs/denmat-energy-acc-v2-20260531-2048-1r
runs/denmat-energy-acc-v2-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 2048 | 1 | 39.00 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy` | 2048 | 1 | 39.40 s | 9.7620 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.88 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy` | 512 | 4 | 9.75 s | 2.1523 GB | 0.000000401 Ha |

| Profile row | 2048/1 baseline | 2048/1 DENMAT GPU | 512/4 baseline, rank 1 | 512/4 DENMAT GPU, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_ETOT_DENMAT` | 3.9086 s | 2.7775 s | 0.2365 s | 0.2003 s |
| `PAW_DENMAT_SITE_KERNEL` | 3.2040 s | 2.0633 s | 0.0816 s | 0.0460 s |
| `PAW_DENMAT_ENERGY_LOOP` | 2.8304 s | 1.6985 s | 0.0792 s | 0.0435 s |
| `ACC_KERNEL_DENMAT_ENERGY_TINV` | - | 0.1265 s | - | 0.0245 s |
| `ACC_COPY_DENMAT_ENERGY_TINV` | - | 4.8633 GB | - | 0.1060 GB |

All runs are energy-valid. This is the first DENMAT GPU path that reduces the
measured DENMAT envelope itself, but the full Si64 wall time is still neutral on
Spark because the prototype transfers Lambda once per site. The next useful
step is therefore not to enable this by default, but to make Lambda/LAGR or the
whole DENMAT working set resident across the atom loop, or to reformulate the
contraction into a small batched BLAS path.

## DENMAT LAGR Residency

The next DENMAT follow-up precomputes `LAGR=LAMBDA*OCC` once per k-point/spin in
the outer `WAVES$DENMAT` loop instead of rebuilding it inside every atom-site
`WAVES_DENMAT` call. When `CPPAW_GPU_DENMAT_ENERGY=1` is active, the same LAGR
block is copied to the GPU once and kept resident across the atom loop. New
profile rows are `PAW_DENMAT_LAGR_SETUP`, `ACC_COPY_DENMAT_LAGR_IN`, and
`ACC_PRESENT_DENMAT_LAGR`.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/denmat-lagr-residency-20260531-512-1r
runs/denmat-lagr-residency-20260531-2048-1r
runs/denmat-lagr-residency-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Total copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 1 | 7.10 s | 1.3676 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy` | 512 | 1 | 6.43 s | 1.3786 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 37.19 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy` | 2048 | 1 | 37.40 s | 4.9892 GB | 0.000000407 Ha |
| `gpu_resident_hpsi` | 512 | 4 | 9.81 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy` | 512 | 4 | 9.56 s | 1.7591 GB | 0.000000401 Ha |

| Profile row | 2048/1 baseline | 2048/1 DENMAT GPU | 512/4 baseline, rank 1 | 512/4 DENMAT GPU, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_DENMAT_LAGR_SETUP` | 0.0217 s | 0.0249 s | 0.0024 s | 0.0022 s |
| `ACC_COPY_DENMAT_LAGR_IN` | - | 0.0758 GB | - | 0.0066 GB |
| `ACC_PRESENT_DENMAT_LAGR` | - | 64 calls | - | 16 calls |
| `PAW_ETOT_DENMAT` | 2.0660 s | 0.8893 s | 0.2280 s | 0.1910 s |
| `PAW_DENMAT_SITE_KERNEL` | 1.3356 s | 0.1399 s | 0.0705 s | 0.0325 s |
| `PAW_DENMAT_ENERGY_LOOP` | 1.3226 s | 0.1268 s | 0.0681 s | 0.0302 s |
| `ACC_KERNEL_DENMAT_ENERGY_TINV` | - | 0.1252 s | - | 0.0214 s |
| `ACC_COPY_DENMAT_ENERGY_TINV` | - | 0.0147 GB | - | 0.0011 GB |
| `PAW_OFFDEN_SUM` | 0.6993 s | 0.7135 s | 0.1524 s | 0.1527 s |

Compared with the previous DENMAT GPU diagnostic, the total copy estimate for
`gpu_resident_hpsi_denmat_energy` drops from 9.7620 GB to 4.9892 GB at 2048/1
and from 2.1523 GB to 1.7591 GB at 512/4. The DENMAT envelope itself also
shrinks strongly, but the full Si64 wall time remains dominated by other phases
at 2048/1. Keep the DENMAT energy path opt-in and use these rows to decide
whether the next step should be a batched BLAS formulation or broader
projector/wavefunction residency around the one-center work.

## Off-Site DENMAT Split

The next instrumentation pass splits the previous `PAW_OFFDEN_SUM` envelope
inside `WAVES_SUMMUPOFFSITEDENMAT` into setup, zeroing, local contraction, and
MPI combine rows:

```
PAW_OFFDEN_SUM_SETUP
PAW_OFFDEN_SUM_ZERO
PAW_OFFDEN_SUM_LOCAL
PAW_OFFDEN_SUM_COMBINE
```

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-profile-split-20260531-512-4r
runs/offden-profile-split-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 4 | 9.75 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy` | 512 | 4 | 9.77 s | 1.7591 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 38.70 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy` | 2048 | 1 | 36.39 s | 4.9892 GB | 0.000000407 Ha |

| Profile row | 2048/1 baseline | 2048/1 DENMAT GPU | 512/4 baseline, rank 1 | 512/4 DENMAT GPU, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_SUM_SETUP` | 0.0000 s | 0.0000 s | 0.0001 s | 0.0000 s |
| `PAW_OFFDEN_SUM_ZERO` | 0.0000 s | 0.0000 s | 0.0001 s | 0.0001 s |
| `PAW_OFFDEN_SUM_LOCAL` | 0.7136 s | 0.6979 s | 0.1386 s | 0.1374 s |
| `PAW_OFFDEN_SUM_COMBINE` | 0.0000 s | 0.0000 s | 0.0152 s | 0.0394 s |

The off-site remainder is therefore mostly local contraction work. At 2048/1
the MPI combine row is negligible; at 512/4 it is visible but still smaller
than the local loop on rank 1. This points the next off-site acceleration pass
toward a local matrix-kernel rewrite or GPU residency around `THIS%PROJ`, not
first toward MPI reduction tuning.

## Off-Site DENMAT BLAS Diagnostic

The follow-up diagnostic adds `CPPAW_GPU_OFFDEN_LOCAL=1` with
`CPPAW_OFFDEN_BLAS=1` as an alias. For scalar `TINV`/`NDIM=1` off-site DENMAT
neighbors it packs the projector factors and evaluates the local contraction
with `ZGEMM`, recording:

```
PAW_OFFDEN_BLAS_PACK
ZGEMM_OFFDEN_TINV_NDIM1
PAW_OFFDEN_BLAS_ACCUM
```

This is still a host-data BLAS prototype, not a fully resident cuBLAS path.
It is intentionally disabled by default. Its role is to measure whether the
local loop has enough matrix-kernel structure to justify a later resident
projector/GPU implementation.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-blas-diagnostic-20260531-512-4r
runs/offden-blas-diagnostic-20260531-2048-1r
runs/offden-blas-combined-20260531-512-4r
runs/offden-blas-combined-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 4 | 9.86 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_blas` | 512 | 4 | 10.16 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi` | 2048 | 1 | 38.79 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_blas` | 2048 | 1 | 37.91 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy` | 512 | 4 | 9.83 s | 1.7591 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 512 | 4 | 9.76 s | 1.7591 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy` | 2048 | 1 | 37.75 s | 4.9892 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 2048 | 1 | 36.91 s | 4.9892 GB | 0.000000407 Ha |

| Profile row | 2048/1 HPSI | 2048/1 HPSI+BLAS | 512/4 HPSI, rank 1 | 512/4 HPSI+BLAS, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_SUM_LOCAL` | 0.7183 s | 0.0723 s | 0.1396 s | 0.0194 s |
| `PAW_OFFDEN_SUM_COMBINE` | 0.0000 s | 0.0000 s | 0.0166 s | 0.0209 s |
| `PAW_OFFDEN_BLAS_PACK` | - | 0.0219 s | - | 0.0058 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | - | 0.0490 s | - | 0.0126 s |
| `PAW_OFFDEN_BLAS_ACCUM` | - | 0.0002 s | - | 0.0001 s |

| Profile row | 2048/1 DENMAT GPU | 2048/1 DENMAT GPU+BLAS | 512/4 DENMAT GPU, rank 1 | 512/4 DENMAT GPU+BLAS, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_SUM_LOCAL` | 0.7151 s | 0.0734 s | 0.1403 s | 0.0181 s |
| `PAW_OFFDEN_SUM_COMBINE` | 0.0000 s | 0.0000 s | 0.0152 s | 0.0158 s |
| `PAW_OFFDEN_BLAS_PACK` | - | 0.0220 s | - | 0.0056 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | - | 0.0498 s | - | 0.0119 s |
| `PAW_OFFDEN_BLAS_ACCUM` | - | 0.0002 s | - | 0.0001 s |

The local contraction shrinks by roughly one order of magnitude in both the
plain HPSI and combined DENMAT-GPU cases, while all energy checks remain within
the existing Si64 tolerance. The total wall-time change is modest because the
Si64 step is dominated by other phases and the prototype still packs host
buffers every neighbor. The next useful GPU step is therefore not another
off-site timing split, but keeping projector/packed off-site buffers resident
and moving this matrix formulation behind an OpenACC/cuBLAS-aware backend.

## Off-Site DENMAT cuBLAS Diagnostic

The follow-up cuBLAS diagnostic adds `CPPAW_GPU_OFFDEN_CUBLAS=1` with
`CPPAW_CUBLAS_ACC_OFFDEN=1` as an alias and a separate
`CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP` threshold. The harness forces that threshold
to 1 for `gpu_resident_*_offden_cublas` cases. This path still packs host
buffers per neighbor, then sends one small `ZGEMM(N,T)` to cuBLAS per neighbor.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-cublas-diagnostic-20260531-512-4r
runs/offden-cublas-diagnostic-20260531-2048-1r
runs/offden-cublas-repro-20260531-2048-1r
runs/offden-cublas-combined-repro-20260531-2048-1r
runs/offden-cublas-combined-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_offden_blas` | 512 | 4 | 10.23 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas` | 512 | 4 | 10.00 s | 1.8762 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_blas` | 2048 | 1 | 40.52 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas` | 2048 | 1 | 39.04 s | 5.3941 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 512 | 4 | 9.75 s | 1.7591 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas` | 512 | 4 | 9.80 s | 1.9069 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 2048 | 1 | 36.79 s | 4.9892 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas` | 2048 | 1 | 37.52 s | 5.4846 GB | 0.000000407 Ha |

| Profile row | 2048/1 HPSI+BLAS | 2048/1 HPSI+cuBLAS | 512/4 HPSI+BLAS, rank 1 | 512/4 HPSI+cuBLAS, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BLAS_PACK` | 0.0227 s | 0.0285 s | 0.0059 s | 0.0083 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | 0.0492 s | - | 0.0118 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_NDIM1` | - | 0.1593 s | - | 0.3990 s |
| `ACC_COPY_CUBLAS_ZGEMM_NT` | - | 0.4954 GB | - | 0.0369 GB |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0733 s | 0.1905 s | 0.0183 s | 0.4084 s |

| Profile row | 2048/1 DENMAT GPU+BLAS | 2048/1 DENMAT GPU+cuBLAS | 512/4 DENMAT GPU+BLAS, rank 1 | 512/4 DENMAT GPU+cuBLAS, rank 1 |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BLAS_PACK` | 0.0222 s | 0.0260 s | 0.0057 s | 0.0080 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | 0.0490 s | - | 0.0118 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_NDIM1` | - | 0.1595 s | - | 0.3998 s |
| `ACC_COPY_CUBLAS_ZGEMM_NT` | - | 0.4954 GB | - | 0.0369 GB |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0727 s | 0.1875 s | 0.0181 s | 0.4090 s |

The naive cuBLAS implementation is therefore a negative control: it is
correct, but the per-neighbor launch/copy pattern loses to host BLAS inside the
off-site local contraction. The useful GPU direction is a batched/strided
formulation grouped by projector shape, ideally with packed `A`/`B` buffers
resident across neighbors. A simple one-cuBLAS-call-per-neighbor rewrite should
not become the default.

## Off-Site DENMAT Stacked cuBLAS Diagnostic

The next diagnostic keeps the same opt-in controls but adds
`CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1` with `CPPAW_CUBLAS_ACC_OFFDEN_BATCH=1` as an
alias. For scalar `TINV`/`NDIM=1`, neighbors with the same first atom and
second-projector size are packed into one wider `ZGEMM(N,T)` call instead of
one tiny cuBLAS call per neighbor. `CPPAW_GPU_OFFDEN_BATCH_SIZE` defaults to 64.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-cublas-stack-diagnostic-20260531-512-4r-fixed
runs/offden-cublas-stack-diagnostic-20260531-2048-1r
runs/offden-cublas-stack-combined-20260531-512-4r
runs/offden-cublas-stack-combined-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_offden_blas` | 512 | 4 | 10.50 s | 1.7284 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas` | 512 | 4 | 10.31 s | 1.8762 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas_batch` | 512 | 4 | 9.65 s | 1.8208 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_blas` | 2048 | 1 | 37.28 s | 4.8988 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas` | 2048 | 1 | 37.47 s | 5.3941 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas_batch` | 2048 | 1 | 37.54 s | 5.1624 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 512 | 4 | 9.78 s | 1.7591 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` | 512 | 4 | 9.84 s | 1.8515 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | 2048 | 1 | 36.31 s | 4.9892 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` | 2048 | 1 | 36.10 s | 5.2528 GB | 0.000000407 Ha |

| Profile row | 2048/1 HPSI+BLAS | 2048/1 HPSI+stacked cuBLAS | 512/4 HPSI+BLAS | 512/4 HPSI+stacked cuBLAS |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BLAS_PACK` | 0.0220 s | - | 0.0247 s | - |
| `PAW_OFFDEN_BATCH_PACK` | - | 0.0692 s | - | 0.0222 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | 0.0486 s | - | 0.0474 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_STACK` | - | 0.0162 s | - | 0.1737 s |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0723 s | 0.0857 s | 0.0746 s | 0.1980 s |

| Profile row | 2048/1 DENMAT GPU+BLAS | 2048/1 DENMAT GPU+stacked cuBLAS | 512/4 DENMAT GPU+BLAS | 512/4 DENMAT GPU+stacked cuBLAS |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BLAS_PACK` | 0.0218 s | - | 0.0235 s | - |
| `PAW_OFFDEN_BATCH_PACK` | - | 0.0655 s | - | 0.0221 s |
| `ZGEMM_OFFDEN_TINV_NDIM1` | 0.0489 s | - | 0.0473 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_STACK` | - | 0.0163 s | - | 0.1739 s |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0727 s | 0.0821 s | 0.0739 s | 0.1981 s |

The stacked form fixes the per-neighbor cuBLAS launch problem: at 2048/1 the
actual GEMM row drops from about 0.157 s in the naive cuBLAS diagnostic to
about 0.016 s. The remaining cost is now host-side packing and copy setup, so
this should stay opt-in. The next implementation step is a resident
projector/off-site buffer path, not forcing this diagnostic as the default.

## Off-Site DENMAT Device-Pack Diagnostic

The follow-up diagnostic adds `CPPAW_GPU_OFFDEN_DEVICE_PACK=1` with
`CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_PACK=1` as an alias. It keeps the stacked
cuBLAS formulation, copies `PROJ` and `OCC` once for the off-site pass, packs
the A/B buffers in OpenACC kernels, calls cuBLAS on present data, and copies
only the stacked `WORK` block back for the existing host-side accumulation.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-device-pack-20260531-512-4r
runs/offden-device-pack-20260531-2048-1r
runs/offden-device-pack-combined-20260531-512-4r
runs/offden-device-pack-combined-20260531-2048-1r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_offden_cublas_batch` | 512 | 4 | 10.60 s | 1.8208 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack` | 512 | 4 | 9.51 s | 1.7485 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas_batch` | 2048 | 1 | 37.45 s | 5.1624 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack` | 2048 | 1 | 37.37 s | 4.9162 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` | 512 | 4 | 10.01 s | 1.8515 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 512 | 4 | 9.87 s | 1.7791 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` | 2048 | 1 | 37.58 s | 5.2528 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 2048 | 1 | 36.03 s | 5.0066 GB | 0.000000407 Ha |

| Profile row | 2048/1 HPSI stacked | 2048/1 HPSI device-pack | 512/4 HPSI stacked | 512/4 HPSI device-pack |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BATCH_PACK` | 0.0692 s | - | 0.0212 s | - |
| `PAW_OFFDEN_DEVICE_PACK` | - | 0.0035 s | - | 0.2960 s |
| `CUBLAS_ZGEMM_OFFDEN_TINV_STACK` | 0.0162 s | - | 0.1738 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_DPACK` | - | 0.0090 s | - | 0.1377 s |
| `ACC_COPY_CUBLAS_ZGEMM_NT` | 0.2636 GB | - | 0.0924 GB | - |
| `ACC_COPY_OFFDEN_DPACK_PROJ_IN` | - | 0.0145 GB | - | 0.0171 GB |
| `ACC_COPY_OFFDEN_DPACK_WORK_OUT` | - | 0.0029 GB | - | 0.0029 GB |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0858 s | 0.0138 s | 0.1969 s | 0.4444 s |

| Profile row | 2048/1 DENMAT GPU stacked | 2048/1 DENMAT GPU device-pack | 512/4 DENMAT GPU stacked | 512/4 DENMAT GPU device-pack |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_BATCH_PACK` | 0.0685 s | - | 0.0225 s | - |
| `PAW_OFFDEN_DEVICE_PACK` | - | 0.0034 s | - | 0.2978 s |
| `CUBLAS_ZGEMM_OFFDEN_TINV_STACK` | 0.0165 s | - | 0.1738 s | - |
| `CUBLAS_ZGEMM_OFFDEN_TINV_DPACK` | - | 0.0090 s | - | 0.1384 s |
| `ACC_COPY_CUBLAS_ZGEMM_NT` | 0.2636 GB | - | 0.0924 GB | - |
| `ACC_COPY_OFFDEN_DPACK_PROJ_IN` | - | 0.0145 GB | - | 0.0171 GB |
| `ACC_COPY_OFFDEN_DPACK_WORK_OUT` | - | 0.0029 GB | - | 0.0029 GB |
| `PAW_OFFDEN_SUM_LOCAL` | 0.0954 s | 0.0137 s | 0.1987 s | 0.4480 s |

Device-pack is the first off-site DENMAT path that materially reduces the
local contraction for the 1 MPI rank / 1 GPU resource comparison. It should not
be promoted to a default for GPU-sharing runs: the 512/4 profile shows the GPU
packing kernels serialize poorly when four ranks share the same device. The
follow-up `THIS%PROJ` residency diagnostic is recorded below. After that, the
next useful implementation step is a device-side accumulation path that avoids
copying `WORK` back for host accumulation.

## `THIS%PROJ` Residency Diagnostic

The follow-up diagnostic adds `CPPAW_GPU_PROJ_RESIDENCY=1` with
`CPPAW_CUBLAS_ACC_PROJ_RESIDENCY=1` as an alias. It copies the combined
`THIS%PROJ` projection result to the device after `WAVES$PROJECTIONS` and lets
eligible off-site DENMAT device-pack paths consume it with `PRESENT_OR_COPYIN`.
Projection recomputation now invalidates any present `PROPSI` output mapping
before entering its `COPYOUT` region; without that guard, later projection calls
could reuse stale device state for an `INTENT(OUT)` result.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/proj-residency-fixed-20260531-2048-1r
runs/proj-residency-fixed-20260531-512-4r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 2048 | 1 | 38.68 s | 5.3749 GB | 0.000000407 Ha |
| `gpu_resident_proj` | 2048 | 1 | 39.44 s | 5.3894 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack` | 2048 | 1 | 38.93 s | 4.9162 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack_proj` | 2048 | 1 | 38.30 s | 4.9162 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 2048 | 1 | 36.93 s | 5.0066 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` | 2048 | 1 | 36.57 s | 5.0066 GB | 0.000000407 Ha |
| `gpu_resident` | 512 | 4 | 9.20 s | 1.8694 GB | 0.000000401 Ha |
| `gpu_resident_proj` | 512 | 4 | 9.27 s | 1.8865 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack` | 512 | 4 | 9.58 s | 1.7485 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_offden_cublas_devicepack_proj` | 512 | 4 | 10.10 s | 1.7484 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 512 | 4 | 9.64 s | 1.7791 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` | 512 | 4 | 9.48 s | 1.7791 GB | 0.000000401 Ha |

| Profile row | 2048/1 no PROJ residency | 2048/1 PROJ residency | 512/4 no PROJ residency | 512/4 PROJ residency |
| --- | ---: | ---: | ---: | ---: |
| `ACC_COPY_THIS_PROJ_IN` | - | 1 call, 0.0145 GB | - | 1 call/rank, 0.0043 GB/rank |
| `ACC_COPY_OFFDEN_DPACK_PROJ_IN` | 1 call, 0.0145 GB | - | 1 call/rank, 0.0043 GB/rank | - |
| `ACC_PRESENT_OFFDEN_DPACK_PROJ` | - | 1 call, 0 GB | - | 1 call/rank, 0 GB |
| `ACC_COPY_PROJ_PROPSI_OUT` | 6 calls, 0.0869 GB | 6 calls, 0.0869 GB | 6 calls/rank, 0.0256 GB/rank | 6 calls/rank, 0.0256 GB/rank |

Conclusion: this is a correctness-safe diagnostic hook for testing broader
projection residency, and it proves the off-site device-pack consumer can avoid
its own `PROJ` copy when the combined projection result is already present. The
net copy estimate does not drop in the current shape because the copy is moved
from the off-site consumer to the post-projection residency step. Keep it
opt-in until more consumers can reuse the same resident `THIS%PROJ` block or the
off-site accumulation path stays fully on the GPU.

## Off-Site DENMAT Flat Device-Accum Diagnostic

The next diagnostic adds `CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1` with
`CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_ACCUM=1` as an alias. It implies the existing
device-pack path, leaves the stacked cuBLAS `ZGEMM` result in device `WORK`,
converts/accumulates each batch into one flat real off-site matrix buffer on the
GPU, and copies that flat buffer back once per k-point/spin pass for the
existing host-side `OSDENMAT` update. This
removes `ACC_COPY_OFFDEN_DPACK_WORK_OUT` and replaces it with
`PAW_OFFDEN_DEVICE_ACCUM`, `ACC_COPY_OFFDEN_DPACK_FLAT_OUT`, and
`PAW_OFFDEN_FLAT_ACCUM_SCATTER`.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/offden-device-accum-20260601-2048-1r
runs/offden-device-accum-20260601-512-4r
runs/offden-flat-accum-20260601-2048-1r
runs/offden-flat-accum-20260601-512-4r
```

| Case | Empty bands | Ranks | Wall time | Copy estimate | Energy delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 2048 | 1 | 37.91 s | 5.0066 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum` | 2048 | 1 | 37.30 s | 5.0052 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` | 2048 | 1 | 37.29 s | 5.0066 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 2048 | 1 | 35.64 s | 5.0052 GB | 0.000000407 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | 512 | 4 | 9.89 s | 1.7791 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum` | 512 | 4 | 9.62 s | 1.7776 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` | 512 | 4 | 9.71 s | 1.7791 GB | 0.000000401 Ha |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 512 | 4 | 9.73 s | 1.7776 GB | 0.000000401 Ha |

| Profile row | 2048/1 device-pack | 2048/1 device-accum | 2048/1 proj device-pack | 2048/1 proj device-accum |
| --- | ---: | ---: | ---: | ---: |
| `PAW_OFFDEN_DEVICE_PACK` | 0.0035 s | 0.0037 s | 0.0036 s | 0.0035 s |
| `CUBLAS_ZGEMM_OFFDEN_TINV_DPACK` | 0.0090 s | 0.0091 s | 0.0091 s | 0.0091 s |
| `ACC_COPY_OFFDEN_DPACK_WORK_OUT` | 0.0029 GB | - | 0.0029 GB | - |
| `PAW_OFFDEN_DEVICE_ACCUM` | - | 0.0006 s, 0.0015 GB | - | 0.0006 s, 0.0015 GB |
| `ACC_COPY_OFFDEN_DPACK_FLAT_OUT` | - | 0.0015 GB | - | 0.0015 GB |
| `PAW_OFFDEN_FLAT_ACCUM_SCATTER` | - | 0.0001 s | - | 0.0001 s |
| `ACC_COPY_OFFDEN_DPACK_PROJ_IN` | 0.0145 GB | 0.0145 GB | - | - |
| `ACC_PRESENT_OFFDEN_DPACK_PROJ` | - | - | 1 call, 0 GB | 1 call, 0 GB |

Conclusion: the implementation is correct and now avoids per-batch result
copies by copying one flat real off-site buffer back at the end of the pass. It
halves the final off-site result copy (`WORK` complex to real flat matrix) from
0.0029 GB to 0.0015 GB and the combined 2048/1 projection-residency case drops
from 37.29 s to 35.64 s. The 512/4 shared-GPU case is neutral/noisy, so this
still belongs behind an opt-in switch. The next useful step is to keep the
consumer of `OSDENMAT` closer to this packed/device representation so the
host-side scatter is no longer the synchronization boundary.

## PSIM Lifecycle Harness

The follow-up adds `tests/profile/si64/run_psim_lifecycle.sh` and
`tests/profile/si64/profile_copy_rows.py`. The lifecycle harness defaults to
`NSTEPS=2`, disables the fixed one-step Si64 energy check, and compares the
reported energies between cases. The copy-row helper reads the profile CSV
header and aggregates the `gbyte` column for `ACC_COPY*` rows, avoiding fragile
positional parsing.

Spark C86C validation:

```
runs/psim-lifecycle-20260601-rerun-empty512-nstep2-1r
runs/psim-lifecycle-20260601-rerun-combined.tsv
runs/psim-lifecycle-20260601-rerun-copy-rows.md
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 2 | 1 | 10.14 s | 2.6165 GB | 269.022536 Ha |
| `gpu_psim_propagate` | 512 | 2 | 1 | 10.05 s | 3.1551 GB | 269.022536 Ha |
| `gpu_resident_psim_phase` | 512 | 2 | 1 | 10.98 s | 2.8861 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 10.05 s | 2.3443 GB | 269.022536 Ha |
| `gpu_hpsi_psim_propagate` | 512 | 2 | 1 | 10.76 s | 2.8829 GB | 269.022536 Ha |
| `gpu_resident_hpsi_psim_phase` | 512 | 2 | 1 | 10.60 s | 2.6139 GB | 269.022536 Ha |

| Profile row | Plain PSIM | PSIM phase | HPSI | HPSI + PSIM phase | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `ACC_COPY_PROP_PSIM_OUT` | 0.1345 GB | - | - | - | Cross-phase PSIM residency removes the immediate propagation copy-out. |
| `ACC_COPY_ORTHO_PSIM_IN` | - | - | 0.1345 GB | - | Plain HPSI still enters orthogonalization through a copied `PSIM`. |
| `ACC_COPY_FORCE_PSI0_IN` | 0.1345 GB | 0.1345 GB | 0.1345 GB | 0.1345 GB | The force path still creates a new `PSI0` residency boundary each step. |
| `ACC_COPY_HPSI_ADDPRO_IN` | - | - | 0.1345 GB | 0.1345 GB | HPSI residency is local to the current ADDPRO/overlap envelope. |
| `ACC_COPY_ADDPRO_OPSI_PSI_IN` | - | 0.1345 GB | 0.1345 GB | 0.1345 GB | The OPSI/ADDPRO side remains a repeated wavefunction copy target. |

Conclusion: the lifecycle run is correctness-clean, but it reinforces the
earlier design choice. PSIM phase residency is doing the intended copy
substitution, yet wall time stays noisy and total copy volume still grows when
PSIM propagation is enabled. The next useful technical lever is broader
`PSI0`/`HPSI`/ADDPRO consumer residency around `WAVES$ETOT` and the following
step, not a standalone promotion of `CPPAW_GPU_PSIM_PHASE_RESIDENCY`.

## OPSI Present-Or-Copy Consumers

The follow-up changes downstream OpenACC data regions that consume already
resident wavefunction buffers from unconditional `copy/copyin` to
`present_or_copy/present_or_copyin`. This covers the fallback projection input
path, the non-cache `WAVES_ADDPRO` input/output region, and both standard and
TINV `WAVES_ADDOPSI` consumers. It deliberately does not keep superwave OPSI
resident through `WAVES_OPSI` build/mass scaling yet, because the superwave
projection path can still read host data.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/opsi-present-consumers-20260601-512-nstep2-1r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 10.22 s | 2.3443 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 10.13 s | 2.3443 GB | 269.022536 Ha |

| Profile row | HPSI | HPSI + OPSI | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_ORTHO_OPSI_IN` | 0.1345 GB | - | Existing OPSI residency removes the later orthogonalization copy-in. |
| `ACC_COPY_OPSI_BUILD_IN` | - | 0.1345 GB | Superwave OPSI still enters residency only after host build/mass scaling. |
| `ACC_COPY_ADDPRO_OPSI_PSI_IN` | 0.1345 GB | 0.1345 GB | The remaining build-time ADDPRO input copy is not solved by consumer cleanup. |
| `ACC_COPY_ADDPRO_OPSI_PSI_OUT` | 0.1345 GB | 0.1345 GB | The build-time ADDPRO output copy remains the next target. |
| `ACC_PRESENT_ADDOPSI_OPSI_TINV` | 2 calls | 2 calls | The TINV ADDOPSI consumer already sees OPSI as present. |
| `ACC_PRESENT_PROJ_ORTHO_OPSI_PSI` | 2 calls | 2 calls | The projection consumer is prepared for resident inputs. |

Conclusion: this is a low-risk enabling patch, not a performance promotion by
itself. It preserves correctness and makes downstream consumers tolerant of
resident inputs, while the measured copy rows point to the next real step:
make the superwave `WAVES_OPSI` build/mass-scale path device-resident only once
the projection path can consume that resident result without falling back to
stale host data.

## Superwave OPSI Build Residency

The follow-up extends `CPPAW_GPU_OPSI_RESIDENCY=1` to the inversion-symmetric
superwave path. `WAVES_OPSI` now starts with OPSI present, `WAVES_ADDPRO`
updates that resident buffer, and mass scaling is performed through the
explicit `WAVES_SCALE_OPSI_ACC` helper. The path still takes one post-mass host
snapshot because later superwave projection fallback code can read OPSI on the
host.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/opsi-superwave-build-residency-20260601-final-v2
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 11.11 s | 2.3443 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 10.36 s | 2.2098 GB | 269.022536 Ha |

| Profile row | HPSI | HPSI + OPSI | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_ADDPRO_OPSI_PSI_IN` | 0.1345 GB | - | OPSI build-time ADDPRO no longer copies the input wavefunction. |
| `ACC_COPY_ADDPRO_OPSI_PSI_OUT` | 0.1345 GB | - | OPSI build-time ADDPRO no longer copies the updated wavefunction back. |
| `ACC_COPY_ORTHO_OPSI_IN` | 0.1345 GB | - | The later orthogonalization consumer sees OPSI present. |
| `ACC_COPY_OPSI_BUILD_IN` | - | 0.1345 GB | OPSI enters the build path on device. |
| `ACC_COPY_OPSI_MASS_OUT` | - | 0.1345 GB | One post-mass host snapshot keeps fallback host readers correct. |
| `ACC_PRESENT_ADDPRO_OPSI_PSI` | - | 2 calls | Build-time ADDPRO updates resident OPSI. |
| `ACC_PRESENT_ORTHO_OPSI` | - | 2 calls | Orthogonalization uses the resident OPSI buffer. |

Conclusion: correctness is preserved for the two-step Si64 smoke and the OPSI
switch now removes the repeated ADDPRO input/output copies for superwave
builds. The remaining host snapshot is intentional; the next cleanup target is
the host-side superwave projection fallback so `ACC_COPY_OPSI_MASS_OUT` can be
reduced or removed later.

## Superwave Projection Residency

The follow-up enables the resident cuBLAS projection path for superwave
wavefunctions. `CPPAW_CUBLAS_ACC_PROJECTION_PRESENT` now applies the same
superwave completion used by `PLANEWAVE$SCALARPRODUCT`: double the projected
sum, subtract the gamma-point term when present, and then apply `GWEIGHT`.
`WAVES_PROJECTIONS` no longer disables the resident path for `SUPER`, and
`CPPAW_GPU_OPSI_RESIDENCY=1` now requires the projection threshold before it
skips the host snapshot after OPSI mass scaling.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/superwave-projection-residency-20260601-512-nstep2-1r
runs/superwave-projection-residency-20260601-512-nstep2-4r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 10.86 s | 2.3443 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 10.61 s | 2.0753 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 4 | 15.24 s | 2.9785 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 4 | 15.14 s | 2.7095 GB | 269.022536 Ha |

| Current OPSI-case profile row | Value | Interpretation |
| --- | ---: | --- |
| `ACC_COPY_OPSI_MASS_OUT` | absent | The superwave projection path no longer needs a post-mass host OPSI snapshot. |
| `ACC_COPY_OPSI_BUILD_IN` | 0.1345 GB | OPSI enters the build path on device. |
| `ACC_PRESENT_ADDPRO_OPSI_PSI` | 2 calls | Build-time ADDPRO updates resident OPSI. |
| `ACC_PRESENT_ORTHO_OPSI` | 2 calls | Orthogonalization sees OPSI present. |
| `ACC_PRESENT_PROJ_ORTHO_OPSI_PSI` | 2 calls | The orthogonalization projection sees OPSI present. |
| `CUBLAS_ZGEMM_PROJ_RES` | 576 calls | Superwave projections use the resident cuBLAS projection kernel. |

Conclusion: the two-step 512-band smoke stays energy-clean for 1 MPI/GPU and
for the 4-rank shared-GPU run. The copy estimate drops by another 0.1345 GB in
the 1-rank case compared with the previous OPSI-build branch, matching removal
of `ACC_COPY_OPSI_MASS_OUT`. This makes the OPSI residency path internally
consistent through build, mass scaling, projection, overlap, and ADDOPSI.

## PSI0/HPSI Copy Boundary

The follow-up keeps the `PSI0` input resident when `CPPAW_GPU_HPSI_RESIDENCY=1`
has already selected the HPSI residency path. The lifetime is deliberately
local to `WAVES$ETOT`: `WAVES$HPSI` creates the device copy if needed, the
immediate expectation and full-Hamiltonian overlap calls reuse it through
OpenACC present checks, and `WAVES$ETOT` deletes that temporary residency before
leaving the energy evaluation.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/psi0-hpsi-copy-boundaries-20260601-512-nstep2-1r
runs/psi0-hpsi-copy-boundaries-20260601-512-nstep2-4r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident` | 512 | 2 | 1 | 11.20 s | 2.6165 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 11.52 s | 2.2098 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 11.76 s | 1.9409 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 4 | 15.46 s | 2.8440 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 4 | 15.33 s | 2.5750 GB | 269.022536 Ha |

| Profile row | HPSI | HPSI + OPSI | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_HPSI_PSI0_IN` | 0.1345 GB | 0.1345 GB | `PSI0` enters the ETOT-local HPSI residency once per step. |
| `ACC_PRESENT_EXPECT_PSI0` | 2 calls | 2 calls | The immediate expectation overlap reuses resident `PSI0`. |
| `ACC_PRESENT_HAMILTON_PSI0` | 2 calls | 2 calls | The full-Hamiltonian overlap reuses resident `PSI0`. |
| `ACC_COPY_EXPECT_PSI0_IN` / `ACC_COPY_HAMILTON_PSI0_IN` | absent | absent | The consumers no longer create their own `PSI0` transfers. |
| `ACC_COPY_HPSI_ADDPRO_IN` | 0.1345 GB | 0.1345 GB | The updated `HPSI` buffer itself still enters the resident addproduct boundary. |
| `ACC_COPY_OPSI_BUILD_IN` | - | 0.1345 GB | OPSI build residency remains explicit and energy-valid. |
| `ACC_PRESENT_ADDPRO_OPSI_PSI` | - | 2 calls | OPSI build-time ADDPRO still updates resident OPSI. |
| `ACC_PRESENT_ORTHO_OPSI` | - | 2 calls | Orthogonalization still sees OPSI present. |

Conclusion: correctness is preserved for both one-rank and four-rank Si64
smokes. Compared with the previous superwave-projection branch, the HPSI and
HPSI+OPSI copy estimates each drop by 0.1345 GB in the one-rank case, exactly
matching the removed duplicate `PSI0` consumer transfer. Wall time remains noisy
at this small size, so this is a copy-boundary cleanup and enabling patch rather
than a new default performance claim.

## Force-To-HPSI PSI0 Residency

The next follow-up reuses the `PSI0` device copy created by the default
force-loop residency when `CPPAW_GPU_HPSI_RESIDENCY=1` is also active.
`WAVES$FORCE` keeps `THIS%PSI0` resident after the per-atom `WAVES_DEDPRO`
loop, and the following `WAVES$HPSI` boundary records that array as present
instead of copying it again. The lifetime still ends inside `WAVES$ETOT` after
the immediate HPSI/expectation/Hamiltonian consumers.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/force-to-hpsi-psi0-residency-20260601-512-nstep2-1r
runs/force-to-hpsi-psi0-residency-20260601-512-nstep2-4r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 11.10 s | 2.0753 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 11.03 s | 1.8064 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 4 | 15.08 s | 2.7095 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 4 | 15.06 s | 2.4405 GB | 269.022536 Ha |

| Profile row | HPSI | HPSI + OPSI | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_FORCE_PSI0_IN` | 0.1345 GB | 0.1345 GB | The force loop creates the `PSI0` device copy once per step. |
| `ACC_PRESENT_HPSI_PSI0` | 2 calls | 2 calls | HPSI reuses the force-loop `PSI0` copy. |
| `ACC_COPY_HPSI_PSI0_IN` | absent | absent | The HPSI-side duplicate `PSI0` copy is removed. |
| `ACC_PRESENT_EXPECT_PSI0` | 2 calls | 2 calls | Expectation still reuses resident `PSI0`. |
| `ACC_PRESENT_HAMILTON_PSI0` | 2 calls | 2 calls | Full-Hamiltonian overlap still reuses resident `PSI0`. |
| `ACC_COPY_HPSI_ADDPRO_IN` | 0.1345 GB | 0.1345 GB | `HPSI` itself still enters the addproduct residency boundary after host-side VPSI/HPROJ work. |
| `ACC_COPY_OPSI_BUILD_IN` | - | 0.1345 GB | OPSI build residency remains the next independent wavefunction copy. |

Conclusion: the copy estimate drops by another 0.1345 GB versus the prior
HPSI/PSI0 boundary branch, with unchanged final energy for 1-rank and 4-rank
smokes. This is still a small-boundary cleanup, but it is a useful step toward
Peter's broader "keep the PAW wavefunctions on the GPU" direction because it
connects two previously separate ETOT-local resident regions without widening
the lifetime beyond the energy evaluation.

## VPSI Device-Finish HPSI Residency

The next boundary cleanup lets `WAVES_VPSI` finish the kinetic/bucket G-space
addition on the device when HPSI residency is active and the input `PSI` is
already present. The FFT/RTOG part is still host-side, so one `HPSI`
host-to-device transfer remains; the point of this patch is to move that copy
accounting to the producing `WAVES_VPSI` boundary and let the following
`WAVES_ADDPRO` step consume `HPSI` as present.
If an HPSI device allocation already exists at this boundary, the refreshed
path records `ACC_UPDATE_VPSI_HPSI_IN` and explicitly updates the device copy
instead of treating the previous device allocation as valid.

The force-to-HPSI `PSI0` carry also uses the same addproduct threshold guard as
the HPSI residency path, so it only extends `PSI0` lifetime when the following
HPSI path can actually use the resident addproduct implementation.

Spark C86C validation:

```
nvhpc_gpu_acc_residency_profile
nvhpc_gpu_acc_residency_profile_parallel

runs/vpsi-device-finish-residency-20260601-512-nstep2-1r
runs/vpsi-device-finish-residency-20260601-512-nstep2-4r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi` | 512 | 2 | 1 | 11.38 s | 2.0755 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 11.11 s | 1.8066 GB | 269.022536 Ha |
| `gpu_resident_hpsi` | 512 | 2 | 4 | 15.16 s | 2.7097 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 4 | 15.36 s | 2.4407 GB | 269.022536 Ha |

| Profile row | HPSI | HPSI + OPSI | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_PRESENT_VPSI_PSI` | 2 calls | 2 calls | VPSI sees the input wavefunction already resident. |
| `ACC_COPY_VPSI_HPSI_IN` | 0.1345 GB | 0.1345 GB | The required HPSI transfer now belongs to the VPSI producer boundary. |
| `ACC_COPY_VPSI_G2_IN` | 0.0002 GB | 0.0002 GB | The kinetic factor copy is tiny compared with the wavefunction transfer. |
| `ACC_PRESENT_HPSI_ADDPRO` | 2 calls | 2 calls | ADDPRO consumes HPSI as present. |
| `ACC_COPY_HPSI_ADDPRO_IN` | absent | absent | The old consumer-side HPSI copy boundary is removed. |
| `ACC_PRESENT_HPSI_PSI0` | 2 calls | 2 calls | HPSI still reuses force-loop `PSI0`. |
| `ACC_PRESENT_EXPECT_PSI0` | 2 calls | 2 calls | Expectation still reuses resident `PSI0`. |
| `ACC_PRESENT_HAMILTON_PSI0` | 2 calls | 2 calls | Full-Hamiltonian overlap still reuses resident `PSI0`. |

Conclusion: this does not yet remove the HPSI host/device copy, because the
FFT/RTOG producer remains host-side. It does make the lifetime graph cleaner:
VPSI produces a device-present `HPSI`, ADDPRO consumes it without another
enter-data boundary, and the profile now points at the real next large target,
namely GPU-resident FFT/RTOG or a broader producer-side wavefunction residency.

## HPSI/OPSI Combined Diagnostic Cases

The follow-up harness patch adds direct `CASES` names for the combinations that
matter now that HPSI and OPSI residency both work:

- `gpu_resident_hpsi_opsi_proj`
- `gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum`
- `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum`

Spark C86C validation used the existing residency binaries:

```
runs/hpsi-opsi-combo-cases-20260601-512-nstep2-1r
runs/hpsi-opsi-combo-cases-20260601-512-nstep2-4r
```

| Case | Empty bands | NSTEPS | Ranks | Wall time | Copy estimate | Final energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 1 | 10.47 s | 1.8066 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_proj` | 512 | 2 | 1 | 11.49 s | 1.8108 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum` | 512 | 2 | 1 | 10.61 s | 1.8181 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 512 | 2 | 1 | 9.58 s | 1.8401 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi` | 512 | 2 | 4 | 15.11 s | 2.4407 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_proj` | 512 | 2 | 4 | 15.03 s | 2.4578 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum` | 512 | 2 | 4 | 15.02 s | 2.4778 GB | 269.022536 Ha |
| `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` | 512 | 2 | 4 | 15.04 s | 2.5391 GB | 269.022536 Ha |

Key profile checks:

| Profile row | Meaning |
| --- | --- |
| `ACC_PRESENT_HPSI_ADDPRO` | HPSI remains present through the ADDPRO consumer. |
| `ACC_PRESENT_ADDPRO_OPSI_PSI` | OPSI build-time ADDPRO updates resident OPSI. |
| `ACC_COPY_THIS_PROJ_IN` plus `ACC_PRESENT_OFFDEN_DPACK_PROJ` | PROJ-residency cases move the projection copy to the producer and let the off-site device-pack consumer see it as present. |
| `PAW_OFFDEN_DEVICE_ACCUM` | Device-pack accumulation remains active when combined with HPSI and OPSI residency. |

Conclusion: the combined switches do not expose a correctness conflict, and the
short Spark smoke is neutral-to-slightly-positive. Treat the full combination
as a diagnostic case for now; the stronger decision still needs the standard
1024/NSTEPS=1 refresh and a larger case before promoting more than HPSI+OPSI
into the recommended path.

## Setup PSI0 Residency Diagnostic

The next setup-boundary check keeps the initial `PSI0` array resident before
`WAVES_GRAMSCHMIDT`, then lets the setup projection and force/HPSI path reuse
that device allocation. The Gram-Schmidt output is still synchronized back to
the host immediately, because later setup and density work still has CPU-side
consumers. The profile therefore counts both the explicit setup upload
(`ACC_COPY_SETUP_PSI0_IN`) and the conservative host update
(`ACC_UPDATE_GRAM_PSI0_PSI_OUT`); this avoids overstating the copy reduction.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries on Spark C86C and Terok:

```
runs/setup-psi0-residency-instr-spark-3rep-20260601-065758
runs/setup-psi0-residency-instr-terok-3rep-20260601-065757
```

| Machine | Case | Empty bands | Repeats | Median wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Spark C86C | `gpu_resident_stack` | 1024 | 3 | 12.18 s | 1.7519 GB | yes |
| Spark C86C | `gpu_resident_stack_setup_host` | 1024 | 3 | 12.35 s | 1.9940 GB | yes |
| Terok A40 | `gpu_resident_stack` | 1024 | 3 | 12.59 s | 1.7519 GB | yes |
| Terok A40 | `gpu_resident_stack_setup_host` | 1024 | 3 | 12.31 s | 1.9940 GB | yes |

Key profile checks from Spark:

| Profile row | New stack | Setup-host control | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_SETUP_PSI0_IN` | 0.3631 GB over 3 repeats | absent | The setup path creates one explicit resident `PSI0` copy. |
| `ACC_UPDATE_GRAM_PSI0_PSI_OUT` | 0.3631 GB over 3 repeats | absent | Gram-Schmidt still refreshes the host copy for CPU consumers. |
| `ACC_COPY_GRAM_PSI0_PSI_IN` | absent | 0.3631 GB over 3 repeats | The Gram input copy is removed. |
| `ACC_COPY_GRAM_PSI0_PSI_OUT` | absent | 0.3631 GB over 3 repeats | The old Gram data-region copyout is replaced by the explicit host update. |
| `ACC_COPY_PROJ_SETUP0_PSI_IN` | absent | 0.3631 GB over 3 repeats | Setup projection reuses resident `PSI0`. |
| `ACC_COPY_FORCE_PSI0_IN` | absent | 0.3631 GB over 3 repeats | The force/HPSI path reuses resident `PSI0`. |

Conclusion: this is correctness-valid and reduces the honest transfer estimate
by about 0.24 GB per 1024-band Si64 step on Spark, but the wall-time signal is
machine/noise dependent. Keep it as part of the residency-stack development
block, not as a standalone performance PR. The next useful optimization is to
shorten or delay the conservative host synchronization only where CPU consumers
can be proven absent, or to attack the still-large `PSIM`, `HPSI`, `OPSI`, and
`WRITEPDOS` wavefunction boundaries.

## Stack Plus PSIM-Phase Diagnostic

The follow-up check combines the focused stack with the existing cross-phase
`PSIM` propagation path via the new harness case
`gpu_resident_stack_psim_phase`. This tests whether `PSIM` should become part of
the stack meta-keyword now that setup `PSI0` residency is in place.

Validation used the already rebuilt residency binaries:

```
runs/stack-psim-phase-spark-1024-nstep1-20260601-071411
runs/stack-psim-phase-terok-1024-nstep1-20260601-071411
```

| Machine | Case | Empty bands | Repeats | Median wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Spark C86C | `gpu_resident_stack` | 1024 | 3 | 12.20 s | 1.7519 GB | yes |
| Spark C86C | `gpu_resident_stack_psim_phase` | 1024 | 3 | 12.13 s | 1.9943 GB | yes |
| Terok A40 | `gpu_resident_stack` | 1024 | 3 | 12.08 s | 1.7519 GB | yes |
| Terok A40 | `gpu_resident_stack_psim_phase` | 1024 | 3 | 12.60 s | 1.9943 GB | yes |

Key copy rows are identical on Spark and Terok for representative repeats:

| Profile row | Stack | Stack + PSIM phase | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_ORTHO_PSIM_IN` | 0.1210 GB | absent | Cross-phase residency removes the orthogonalization input copy. |
| `ACC_COPY_PROP_PSI0_IN` | absent | 0.1210 GB | Propagation now needs a `PSI0` input on device. |
| `ACC_COPY_PROP_PSIM_IN` | absent | 0.1210 GB | Propagation also needs the previous `PSIM`. |
| `ACC_COPY_PROP_HPSI_IN` | absent | 0.1210 GB | Propagation also needs `HPSI`. |
| `ACC_COPY_ORTHO_PSIM_OUT` | 0.1210 GB | 0.1210 GB | The final host refresh remains. |

Conclusion: this is a useful diagnostic and remains energy-valid, but it should
not be folded into `CPPAW_GPU_RESIDENCY_STACK` yet. It trades one removed
orthogonalization input copy for three propagation input copies, increasing the
honest transfer estimate by about 0.24 GB per 1024-band step; Spark timing is
noise-level neutral and Terok is worse in the median.

## PSI0 To OPSI Residency Diagnostic

The next residency step keeps `PSI0` resident across the ETOT cleanup boundary
only when the OPSI residency path is enabled, then builds `OPSI` directly on
the device in `WAVES$ORTHOGONALIZE`. The behavior is controlled by
`CPPAW_GPU_PSI0_ORTHO_RESIDENCY`; it is enabled by the stack meta-keyword and
can be disabled with `CPPAW_GPU_PSI0_ORTHO_RESIDENCY=0`.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries on Spark C86C and Terok:

```
runs/psi0-ortho-residency-spark-r3-20260601-073321
runs/psi0-ortho-residency-terok-r3-20260601-073321
runs/psi0-ortho-control-spark-smoke-20260601-073508
runs/psi0-ortho-control-terok-smoke-20260601-073508
```

| Machine | Case | Empty bands | Repeats | Median wall time | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Spark C86C | `gpu_resident_stack` | 1024 | 3 | 12.18 s | 1.6309 GB | yes |
| Spark C86C | `gpu_resident_stack_psi0_ortho_host` | 1024 | 1 | 12.37 s | 1.7519 GB | yes |
| Terok A40 | `gpu_resident_stack` | 1024 | 3 | 12.16 s | 1.6309 GB | yes |
| Terok A40 | `gpu_resident_stack_psi0_ortho_host` | 1024 | 1 | 12.22 s | 1.7519 GB | yes |

Key profile checks are identical on Spark and Terok for representative
repeats:

| Profile row | New stack | PSI0-ortho-host control | Interpretation |
| --- | ---: | ---: | --- |
| `ACC_COPY_OPSI_BUILD_IN` | absent | 0.1210 GB | The old host-built OPSI copy is removed. |
| `ACC_PRESENT_OPSI_BUILD_PSI0` | present | absent | OPSI now starts from resident `PSI0`. |
| `ACC_CREATE_OPSI_BUILD` | present | absent | OPSI is allocated on the device before the copy kernel. |
| `ACC_COPY_SETUP_PSI0_IN` | 0.1210 GB | 0.1210 GB | Setup still creates the resident `PSI0` allocation. |
| `ACC_UPDATE_GRAM_PSI0_PSI_OUT` | 0.1210 GB | 0.1210 GB | The conservative Gram-Schmidt host refresh remains. |

Conclusion: this is correctness-valid on both GPU machines and removes one
1024-band wavefunction transfer, reducing the honest transfer estimate by
0.1210 GB per Si64 step. Wall time remains noise-level neutral, but this is the
first positive cross-boundary residency result that directly supports the
"keep wavefunction data on the GPU" direction.

## Single-Rank Full-Grid 3D FFT Diagnostic

`CPPAW_FFT_SERIAL_3D=1` adds an opt-in single-rank `PLANEWAVE$FFT` path that
maps the plane-wave stripes to a full local 3D grid, calls `LIB$3DFFTC8`, and
maps the result back. In GPU builds this can be combined with
`CPPAW_CUFFT_ACC=1 CPPAW_CUFFT_ACC_3D=1 CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS=0`;
the harness case is `gpu_resident_stack_serial3dfft`.

Validation used rebuilt `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` binaries on Spark C86C:

```
runs/serial3dfft-probe-20260601-131207
runs/serial3dfft-nsteps3-20260601-131244
runs/serial3dfft-parallel-fallback-20260601-131530
```

| Machine | Case | Empty bands | NSTEPS | Wall time | `vpsi_s` | `fft_s` | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Spark C86C | `gpu_resident_stack` | 1024 | 1 | 12.26 s | 1.3298 s | 3.8647 s | 1.5102 GB | yes |
| Spark C86C | `gpu_resident_stack_serial3dfft` | 1024 | 1 | 10.56 s | 0.3453 s | 0.9354 s | 5.6070 GB | yes |
| Spark C86C | `gpu_resident_stack` | 1024 | 3 | 28.40 s | 3.9662 s | 11.4046 s | 2.9994 GB | check disabled |
| Spark C86C | `gpu_resident_stack_serial3dfft` | 1024 | 3 | 23.53 s | 1.0505 s | 2.6596 s | 15.2897 GB | check disabled |

Representative profile rows for the one-step case:

| Profile row | Stack | Serial 3D FFT | Interpretation |
| --- | ---: | ---: | --- |
| `PAW_VPSI_TOTAL` | 1.3298 s | 0.3453 s | The local 3D cuFFT path removes most of the stripe FFT envelope cost for one rank. |
| `PW_FFT_SERIAL3D_TOTAL` | absent | 0.6539 s | New full-grid path covers the `PLANEWAVE$FFT` calls. |
| `CUFFT3D_C8` | absent | 0.2590 s | Actual cuFFT time is a fraction of the serial-3D envelope. |
| `ACC_COPY_CUFFT3D_C8` | absent | 4.0968 GB | The current wrapper still copies full grids in/out. |

The 4-rank fallback smoke (`EMPTY_BANDS=128`, `NSTEPS=1`) is energy-valid and
shows no `PW_FFT_SERIAL3D_TOTAL`/`CUFFT3D_C8` rows, so the new path is confined
to `NTASKS=1`. Conclusion: this is the first FFT-side GPU diagnostic with a
clear wall-time win for the intended one-rank/one-GPU comparison, despite the
larger explicit transfer estimate. It should remain opt-in until a resident
full-grid cuFFT path removes the copy volume.

## Force DEDPRO/PROFORCE Residency Diagnostic

`CPPAW_GPU_FORCE_DEDPRO_RESIDENCY=1` adds an opt-in force-loop diagnostic for
non-stress `TINV`, `NDIM=1` cases. It keeps the cuBLAS-built `DEDPRO` matrix on
the GPU, applies the inversion-symmetry averaging on device, and evaluates the
`WAVES_PROFORCE` contraction as GPU reductions so only the three force
components return to the host. The harness case is
`gpu_resident_stack_force_dedpro`.

Validation rebuilt both `nvhpc_gpu_acc_residency_profile` and
`nvhpc_gpu_acc_residency_profile_parallel` on Spark C86C and Terok. One-rank
`NSTEPS=1` runs passed the Si64 energy check on both machines, and four-rank
smokes also passed. The longer `NSTEPS=3` diagnostic used the energy check
disabled because the expected final energy differs from the one-step reference.

| Machine | Case | Empty bands | NSTEPS | Wall time | Transfer estimate | Energy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Spark C86C | `gpu_resident_stack` | 1024 | 1 | 12.12 s | 1.3892 GB | 302.280854 |
| Spark C86C | `gpu_resident_stack_force_dedpro` | 1024 | 1 | 12.25 s | 1.3455 GB | 302.280854 |
| Terok | `gpu_resident_stack` | 1024 | 1 | 12.58 s | 1.3892 GB | 302.280853 |
| Terok | `gpu_resident_stack_force_dedpro` | 1024 | 1 | 12.45 s | 1.3455 GB | 302.280853 |
| Spark C86C | `gpu_resident_stack` | 1024 | 3 | 28.08 s | 2.8783 GB | 208.886424 |
| Spark C86C | `gpu_resident_stack_force_dedpro` | 1024 | 3 | 27.46 s | 2.7472 GB | 208.886424 |
| Terok | `gpu_resident_stack` | 1024 | 3 | 28.96 s | 2.8783 GB | 208.886424 |
| Terok | `gpu_resident_stack_force_dedpro` | 1024 | 3 | 28.91 s | 2.7472 GB | 208.886424 |

Representative `NSTEPS=3` profile rows show the intended transfer shift:
`ACC_COPY_ZGEMM_MAT_C_OUT` disappears from the force `DEDPRO` calls, while
`ACC_COPY_FORCE_DEDPRO_INPUTS_IN`, `CUBLAS_ZGEMM_FORCE_DEDPRO`, and
`ACC_FORCE_PROFORCE` appear. The transfer saving is stable at about 44 MB per
Si64 step with 1024 empty bands. Wall time is modestly positive on Spark and
neutral on Terok, so this stays opt-in rather than joining the default
`CPPAW_GPU_RESIDENCY_STACK`.

The Spark follow-up also checked the force path together with the single-rank
3-D cuFFT diagnostic. Runs were sequential on the same GPU; parallel probe runs
were ignored because they contended for the device.

| Case | Empty bands | NSTEPS | Wall time | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident_stack_serial3dfft` | 1024 | 1 | 11.17 s | 5.4859 GB | yes |
| `gpu_resident_stack_serial3dfft_force_dedpro` | 1024 | 1 | 10.63 s | 5.4422 GB | yes |
| `gpu_resident_stack_serial3dfft` | 1024 | 3 | 23.70 s | 15.1687 GB | check disabled |
| `gpu_resident_stack_serial3dfft_force_dedpro` | 1024 | 3 | 23.38 s | 15.0375 GB | check disabled |

Conclusion: the force DEDPRO device path remains a small opt-in diagnostic on
its own, but it composes cleanly with the much stronger one-rank 3-D cuFFT path
and saves the same force-transfer volume there. The harness case is
`gpu_resident_stack_serial3dfft_force_dedpro`.

## Capability-Driven Overnight Smoke

The capability-driven overnight harness was refreshed on 2026-06-01 after the
default case selection was moved to `paw_gpu_capabilities.sh`. Dry-runs on Spark
C86C and Terok both selected the same resource comparison:
`cpu nvhpc_cpu gpu_resident_stack` where a GNU CPU binary exists, plus focused
GPU diagnostics, band cases, and the `gpu_resident_stack` Nsight target. Dry-run
mode now also skips Nsight collection instead of accidentally starting `nsys`.

Spark C86C then ran a real `NSTEPS=1` smoke after building the missing
`nvhpc_gpu_acc_residency_profile_parallel` target. Run root:
`tests/profile/si64/runs/overnight-smoke-spark-20260601-172446`.

| Suite | Case | Ranks | Wall time | Energy delta |
| --- | --- | ---: | ---: | ---: |
| `main_1steps_4ranks` | `cpu` | 4 | 3.75 s | 0.000000 |
| `main_1steps_4ranks` | `nvhpc_cpu` | 4 | 3.84 s | 0.000000 |
| `main_1steps_4ranks` | `gpu_resident_stack` | 4 | 5.25 s | 0.000000 |
| `scaling_1steps_1ranks` | `cpu` | 1 | 3.38 s | 0.000000 |
| `scaling_1steps_1ranks` | `nvhpc_cpu` | 1 | 3.40 s | 0.000000 |
| `scaling_1steps_1ranks` | `gpu_resident_stack` | 1 | 2.95 s | 0.000000 |
| `scaling_1steps_2ranks` | `cpu` | 2 | 5.54 s | 0.000000 |
| `scaling_1steps_2ranks` | `nvhpc_cpu` | 2 | 5.58 s | 0.000000 |
| `scaling_1steps_2ranks` | `gpu_resident_stack` | 2 | 4.81 s | 0.000000 |
| `scaling_1steps_4ranks` | `cpu` | 4 | 3.71 s | 0.000000 |
| `scaling_1steps_4ranks` | `nvhpc_cpu` | 4 | 3.83 s | 0.000000 |
| `scaling_1steps_4ranks` | `gpu_resident_stack` | 4 | 5.14 s | 0.000000 |
| `threshold_1e7_1steps_4ranks` | `gpu_resident_stack` | 4 | 4.98 s | 0.000000 |

Terok exposed a useful installation bug: the first real smoke found the NVHPC
SDK under `$HOME/opt/nvidia/hpc_sdk`, but the runtime launcher search only
checked `$NVHPC_ROOT` and `/opt/nvidia/...`. Serial runs passed, while all MPI
runs failed with `timeout: failed to run command 'mpirun': No such file or
directory`. The harness now searches `$HOME/opt/nvidia/hpc_sdk/<platform>/*`
and adjacent NVHPC version directories, and `paw_gpu_capabilities.sh` reports
the resolved `mpirun` and `nvfortran` paths. The fixed run root:
`tests/profile/si64/runs/overnight-smoke-terok-mpirunfix-20260601-173611`.

| Suite | Case | Ranks | Wall time | Energy delta |
| --- | --- | ---: | ---: | ---: |
| `main_1steps_4ranks` | `nvhpc_cpu` | 4 | 3.10 s | 0.000001 |
| `main_1steps_4ranks` | `gpu_resident_stack` | 4 | 4.91 s | 0.000001 |
| `scaling_1steps_1ranks` | `nvhpc_cpu` | 1 | 4.94 s | 0.000001 |
| `scaling_1steps_1ranks` | `gpu_resident_stack` | 1 | 3.83 s | 0.000001 |
| `scaling_1steps_2ranks` | `nvhpc_cpu` | 2 | 4.10 s | 0.000001 |
| `scaling_1steps_2ranks` | `gpu_resident_stack` | 2 | 4.21 s | 0.000001 |
| `scaling_1steps_4ranks` | `nvhpc_cpu` | 4 | 3.11 s | 0.000001 |
| `scaling_1steps_4ranks` | `gpu_resident_stack` | 4 | 4.63 s | 0.000001 |
| `threshold_1e7_1steps_4ranks` | `gpu_resident_stack` | 4 | 4.57 s | 0.000001 |

Conclusion: the one-rank GPU path is consistently faster than one-rank CPU on
both hosts, while the four-rank CPU/NVHPC path remains faster for Si64. This
keeps the design direction unchanged: Si64 is a correctness and harness smoke,
not the deciding performance target; larger band/projection-heavy cases remain
necessary before promoting more GPU residency paths into defaults.

## 2026-06-01 2048-Band Resource Comparison

The larger follow-up uses the same current branch and binaries, but increases
the Si64 band stress to `EMPTY_BANDS=2048`, keeps `NSTEPS=1`, and compares one
MPI rank plus one GPU against one-rank and eight-rank `nvhpc_cpu`. The GPU list
is focused on the current resident stack and the strongest single-rank FFT/force
diagnostics, without the fallback/off matrix:

```
GPU_CASES="gpu_resident_stack gpu_resident_stack_serial3dfft \
gpu_resident_stack_serial3dfft_force_dedpro \
gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache"
CPU_CASES="nvhpc_cpu" GPU_RANKS=1 CPU_RANKS=8 ADD_GPU_FALLBACK=no
```

Run directories:

```
Spark: tests/profile/si64/runs/si64-bands2048-focused-spark-20260601-174042
Terok: tests/profile/si64/runs/si64-bands2048-focused-terok-20260601-174041
```

| System | Best one-GPU case | `gpu_resident_stack` | 1-rank `nvhpc_cpu` | 8-rank `nvhpc_cpu` | Best GPU speedup | Energy delta |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Spark GB10 | `gpu_resident_stack_serial3dfft_force_dedpro` 32.06 s | 35.28 s | 388.98 s | 1122.63 s | 12.13x vs 1 CPU, 35.02x vs 8 CPU | 0.0000004 |
| Terok A40 | `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache` 32.79 s | 36.69 s | 1109.15 s | 1341.45 s | 33.83x vs 1 CPU, 40.91x vs 8 CPU | 0.0000009 |

GPU case details:

| Case | Spark wall / transfer | Terok wall / transfer | Interpretation |
| --- | ---: | ---: | --- |
| `gpu_resident_stack` | 35.28 s / 2.95 GB | 36.69 s / 2.95 GB | Conservative default remains strong and portable. |
| `gpu_resident_stack_serial3dfft` | 32.40 s / 10.69 GB | 34.06 s / 10.69 GB | Faster on both systems despite the explicit full-grid copy volume. |
| `gpu_resident_stack_serial3dfft_force_dedpro` | 32.06 s / 10.64 GB | 34.89 s / 10.64 GB | Best on Spark; force DEDPRO composes cleanly with serial 3-D FFT. |
| `gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_cache` | 42.92 s / 4.01 GB | 32.79 s / 4.01 GB | Best on Terok, but still hurts Spark; keep ACCMAP/cache opt-in. |

This is the first larger resource comparison where the GPU result is not merely
faster than a serial CPU reference: on both machines, one MPI rank with one GPU
also beats the eight-rank CPU/NVHPC run decisively for this single-k-point,
band-heavy Si64 stress. The default should therefore stay with the conservative
resident stack, while the single-rank serial 3-D FFT and force-DEDPRO paths are
worth keeping as explicit benchmark cases. The ACCMAP/cache path remains a
system-dependent diagnostic until its Spark behavior is understood.

## 2026-06-01 ACCMAP Split/Cache Validation

Commit `b4a18ed` splits the serial 3-D ACCMAP transfer accounting into
`ACC_COPY_SERIAL3D_ACC_INPUT`, `ACC_COPY_SERIAL3D_ACC_OUTPUT`, and
`ACC_COPY_SERIAL3D_ACC_MAP_META`, while keeping the legacy total row
`ACC_COPY_SERIAL3D_ACC_MAP`. It also adds
`gpu_resident_stack_serial3dfft_accmap_cache` so the ACCMAP cache can be tested
without the HPSI/VPSI residency follow-ups.

Run directories:

```
Spark: tests/profile/si64/runs/accmap-split-cache-spark-20260601-185329
Terok: tests/profile/si64/runs/accmap-split-cache-terok-20260601-185329
```

Both systems used `EMPTY_BANDS=2048`, `NSTEPS=1`, one MPI rank and one GPU.

| System | `*_accmap_cache` | `*_hpsi_rtog_vpsi_internal_cache` | Energy delta |
| --- | ---: | ---: | ---: |
| Spark GB10 | 43.13 s / 7.05 GB | 42.35 s / 4.01 GB | 0.000000 |
| Terok A40 | 33.48 s / 7.05 GB | 32.19 s / 4.01 GB | 0.000001 |

The split rows show the same byte accounting on both systems. For the isolated
cache case, the serial-3D transfer total is 4.0978 GB, split into 1.2897 GB
input, 2.8081 GB output, and 0.0000 GB map metadata. For the combined
HPSI/VPSI residency cache case, input drops to 0.0000 GB, output remains
1.2897 GB, map metadata remains 0.0000 GB, and the legacy total row is
1.2897 GB. That confirms the cache removes repeated map metadata traffic and
the combined residency removes the remaining input traffic; the main remaining
ACCMAP transfer is the RTOG/GTOR output boundary. The split rows are treated as
diagnostic detail rows in the summary tools so the legacy total row remains the
only serial-3D ACCMAP contribution to the aggregate transfer estimate.

Interpretation is unchanged but sharper: Terok still benefits from the full
ACCMAP/cache/residency combination, while Spark still pays a runtime cost even
after the byte volume is lower. Keep ACCMAP/cache opt-in and use the split rows
to decide whether the next useful step is output-boundary residency or avoiding
the serial mapping kernels entirely.

## 2026-06-01 ACCMAP Output Present-Or-Copyout Follow-Up

Commit `e475c75` changes the serial 3-D ACCMAP data region from unconditional
`COPYOUT` to `PRESENT_OR_COPYOUT` for the GTOR/RTOG output arrays. This keeps
the fallback semantics when no caller-owned device output exists, but avoids
forcing an output boundary when a resident caller buffer can be used.

Run directories:

```
Spark: tests/profile/si64/runs/accmap-presentout-spark-20260601-190535
Terok: tests/profile/si64/runs/accmap-presentout-terok-20260601-190534
```

Both systems used `EMPTY_BANDS=2048`, `NSTEPS=1`, one MPI rank and one GPU.

| System | `*_accmap_cache` | `*_hpsi_rtog_vpsi_internal_cache` | Energy delta |
| --- | ---: | ---: | ---: |
| Spark GB10 | 41.32 s / 7.05 GB | 42.12 s / 4.01 GB | 0.000000 |
| Terok A40 | 33.58 s / 7.05 GB | 31.89 s / 4.01 GB | 0.000001 |

Compared with the previous split/cache validation, the isolated cache case
improves on Spark from 43.13 s to 41.32 s and is noise-neutral on Terok
(33.48 s to 33.58 s). The full ACCMAP/HPSI/VPSI cache case improves slightly on
both systems: Spark 42.35 s to 42.12 s, Terok 32.19 s to 31.89 s. The split
transfer rows are unchanged because they are conservative caller-boundary
accounting rows; runtime is the relevant evidence here.

This is worth keeping because it is semantically narrower than a new residency
mode and composes cleanly with existing cases, but it does not change the
default recommendation: ACCMAP/cache still remains opt-in until Spark's
mapping-kernel/runtime cost is reduced.

## 2026-06-01 ACCMAP Resident Data-Region Skip

Commit `32b4d4d` adds a resident fast path to the serial 3-D ACCMAP data
region. If the input array, output array, cached full-grid work array, and map
metadata are already present on the device, the structured OpenACC data region
is skipped and the profiler records `ACC_PRESENT_SERIAL3D_DATA`. Non-resident
callers continue to use the same `PRESENT_OR_*` fallback clauses.

Run directories:

```
Spark: tests/profile/si64/runs/accmap-dataskip-spark-20260601-191252
Terok: tests/profile/si64/runs/accmap-dataskip-terok-20260601-191252
```

Both systems used `EMPTY_BANDS=2048`, `NSTEPS=1`, one MPI rank and one GPU.

| System | `*_accmap_cache` | `*_hpsi_rtog_vpsi_internal_cache` | Energy delta |
| --- | ---: | ---: | ---: |
| Spark GB10 | 44.17 s / 7.05 GB | 42.17 s / 4.01 GB | 0.000000 |
| Terok A40 | 33.13 s / 7.05 GB | 31.68 s / 4.01 GB | 0.000001 |

The new present row appears only in the fully resident combined case, with
2176 `ACC_PRESENT_SERIAL3D_DATA` calls on both systems. A follow-up direction
split, validated after removing an experimental caller hint, resolves that
total into 1088 `ACC_PRESENT_SERIAL3D_GTOR_DATA` and 1088
`ACC_PRESENT_SERIAL3D_RTOG_DATA` calls:

```
Spark: tests/profile/si64/runs/accmap-direction-spark-20260601-192642
Terok: tests/profile/si64/runs/accmap-direction-terok-20260601-192641
```

That is the expected shape: the isolated cache case still has a non-resident
input/output boundary, while `*_hpsi_rtog_vpsi_internal_cache` keeps the
`WAVES_VPSI` GTOR/RTOG pair resident enough to skip the data region. The
remaining 1.2897 GB `ACC_COPY_SERIAL3D_ACC_OUTPUT` in the combined case is
therefore not the VPSI RTOG boundary; it comes from another GTOR output
boundary that still lacks a resident consumer.

The change is correctness-safe in this smoke case, but its timing impact is
system-dependent. Terok's `PW_FFT_SERIAL3D_TOTAL` drops from about 0.87 s in the
previous present-or-copyout run to 0.74 s, and the wall time improves slightly
from 31.89 s to 31.68 s. Spark confirms that this particular data-region
boundary is not the remaining large cost: `PW_FFT_SERIAL3D_TOTAL` stays near
9.87 s and wall time remains essentially unchanged for the combined case
(42.12 s to 42.17 s). Keep the skip because it removes avoidable runtime
bookkeeping in the resident path, but keep ACCMAP/cache opt-in and continue
looking for Spark's remaining serial-3D overhead elsewhere.

## 2026-06-01 Density GTOR Residency Prototype

Commit `7b19235` adds the opt-in
`CPPAW_GPU_DENSITY_INTERNAL_RESIDENCY=1` diagnostic. For scalar,
non-kinetic-density calls in `WAVES_DENSITY`, the GTOR output stays resident on
the GPU, `RHO` is accumulated with OpenACC, and only the final density is copied
back. The harness case is
`gpu_resident_stack_serial3dfft_accmap_hpsi_rtog_vpsi_internal_density_cache`.

Run directories:

```
Spark: tests/profile/si64/runs/density-spark-20260601-194229
Terok: tests/profile/si64/runs/density-terok-20260601-194229
```

Both systems used `EMPTY_BANDS=2048`, `NSTEPS=1`, one MPI rank and one GPU.

| System | Cache baseline | + density residency | Energy delta |
| --- | ---: | ---: | ---: |
| Spark GB10 | 40.43 s / 4.01 GB | 31.49 s / 2.72 GB | 0.000000 |
| Terok A40 | 32.00 s / 4.01 GB | 31.54 s / 2.72 GB | 0.000001 |

The transfer reduction is exactly the formerly remaining density-side serial
3-D output boundary: `ACC_COPY_SERIAL3D_ACC_OUTPUT` drops from 1.2897 GB to
zero bytes, while the new diagnostic records one
`ACC_PRESENT_DENS_PSIOFR(74088,1,1088)` row. The only new explicit density
transfers are small: occupations (`0.000017 GB`), `RHO` copyout
(`0.000593 GB`), and the time-inversion phase factor (`0.001185 GB`).

The density timing confirms the source of Spark's large serial-3D overhead. On
Spark, `PAW_DENSITY_TOTAL` drops from 8.9750 s to 0.1588 s
(`PAW_DENSITY_FFT`: 8.9092 s to 0.1555 s). On Terok, the same path drops from
0.4697 s to 0.1044 s. Energy remains within the existing tolerance on both
systems. This is the strongest current evidence that broader caller-side
residency around `PLANEWAVE$FFT` is more valuable than tuning the isolated
ACCMAP kernels alone, especially on Spark.

## 2026-06-01 Batched Orthogonalization 1COVERLAP

The opt-in `CPPAW_GPU_1COVERLAP_BATCH=1` diagnostic adds a combined cuBLAS path
for the three `WAVES_1COVERLAP` calls inside `WAVES$ORTHOGONALIZE`: `PROJ`,
`OPROJ`, `DO*PROJ`, and `DO*OPROJ` are packed once, then the three overlap
matrices are formed in one device data region. The default remains unchanged.

Run directories:

```
Spark: tests/profile/si64/runs/1cov-batch-spark-repeat-20260601-200631
Spark 4-rank smoke: tests/profile/si64/runs/1cov-batch-spark-4r-20260601-201111
Terok: tests/profile/si64/runs/1cov-batch-terok-20260601-201409
```

The large one-rank comparisons used `EMPTY_BANDS=2048`, `NSTEPS=1`.

| System | Case | Wall time | Transfer estimate | Projector copy estimate | Energy check |
| --- | --- | ---: | ---: | ---: | --- |
| Spark GB10 median of 3 | density stack | 32.58 s | 2.7218 GB | 0.3476 GB | yes |
| Spark GB10 median of 3 | density stack + 1C batch | 32.09 s | 2.6494 GB | 0.2752 GB | yes |
| Terok A40 single run | density stack | 31.14 s | 2.7218 GB | 0.3476 GB | yes |
| Terok A40 single run | density stack + 1C batch | 30.93 s | 2.6494 GB | 0.2752 GB | yes |

The parallel smoke used `EMPTY_BANDS=512`, `NSTEPS=1`, and four MPI ranks on
Spark. It compares the focused one-center overlap diagnostic without density
residency, because the density prototype intentionally requires the serial
ACCMAP GTOR path.

| Case | Ranks | Wall time | Transfer estimate | Projector copy estimate | Energy check |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident_1coverlap` | 4 | 9.44 s | 1.8694 GB | 0.2002 GB | yes |
| `gpu_resident_1coverlap_batch` | 4 | 9.31 s | 1.8482 GB | 0.1789 GB | yes |

Representative Spark profile rows show that the orthogonalization part moves
from three legacy `WAVES_1COVERLAP` calls to one batch row, while the two
Gram-Schmidt calls still use the existing helper. In the single-run profile,
`PAW_ORTHO_1COVERLAP` is essentially neutral (`0.2523 s` baseline vs
`0.2451 s` batch), and the real value is cleaner accounting plus a modest
projector-transfer reduction. Keep the switch opt-in for now; it is correct and
slightly favorable in these runs, but the total wall-time gain is well within
run-to-run noise.

## 2026-06-01 Orthogonalization ADDOPROJ Slice GEMMs

The opt-in `CPPAW_GPU_ORTHO_ADDOPROJ=1` diagnostic moves the large
orthogonalization `WAVES_ADDOPROJ` update to cuBLAS slice GEMMs:
`PROJ(:,:,p) += OPROJ(:,:,p) * LAMBDA`. The time-inversion path builds
`LAMBDA1`/`LAMBDA2` on the host, copies them once for the large projector block,
and applies both the direct and conjugated `OPROJ` contributions on the GPU.
`CPPAW_GPU_ORTHO_ADDOPROJ_MIN_NPRO` defaults to `64`; this avoids offloading the
small per-atom `NPRO=13` calls that otherwise copy the same large Lambda blocks
many times.

Run directories:

```
Spark: tests/profile/si64/runs/addoproj-threshold-2048-20260601
Spark 4-rank smoke: tests/profile/si64/runs/addoproj-threshold-4r-20260601
Terok: tests/profile/si64/runs/addoproj-threshold-terok-2048-20260601
```

The large one-rank comparisons used `EMPTY_BANDS=2048`, `NSTEPS=1`, the
density-resident stack, and the batched one-center overlap path.

| System | Case | Wall time | `PAW_ORTHO_ADDOPROJ` | Transfer estimate | Energy check |
| --- | --- | ---: | ---: | ---: | --- |
| Spark GB10 median of 3 | density stack + 1C batch | 33.07 s | 2.6095 s | 2.6494 GB | yes |
| Spark GB10 median of 3 | + ADDOPROJ slice GEMMs | 30.13 s | 0.1610 s | 2.7307 GB | yes |
| Terok A40 single run | density stack + 1C batch | 31.98 s | 2.1271 s | 2.6494 GB | yes |
| Terok A40 single run | + ADDOPROJ slice GEMMs | 29.74 s | 0.1106 s | 2.7307 GB | yes |

The added transfer is small after the `NPRO` threshold: Spark's projector copy
estimate rises from `0.2752 GB` to `0.3565 GB`, mostly from one copy each of
`PROJ`, `OPROJ`, `LAMBDA1`, and `LAMBDA2`. The initial ungated prototype also
offloaded 64 small per-atom calls and raised total transfers to `5.1985 GB`;
the threshold removes that side effect while preserving the large-block speedup.

The four-rank Spark smoke used `EMPTY_BANDS=512`, `NSTEPS=1`, and the focused
non-density cases because density residency is intentionally restricted to the
serial 3-D ACCMAP path.

| Case | Ranks | Wall time | `PAW_ORTHO_ADDOPROJ` | Transfer estimate | Energy check |
| --- | ---: | ---: | ---: | ---: | --- |
| `gpu_resident_1coverlap_batch` | 4 | 9.31 s | 2.0954 s rank-summed | 1.8482 GB | yes |
| `gpu_resident_addoproj` | 4 | 8.94 s | about 0.39 s rank-summed | 1.9337 GB | yes |

This is the first orthogonalization-side diagnostic in the current branch with
a clear wall-time signal beyond copy accounting. Keep it opt-in until a larger
multi-step case confirms that the extra projector/Lambda traffic remains small
outside the Si64 one-step benchmark.

## 2026-06-01 cuSOLVER Gram-Cholesky Diagnostic

`CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY=1` enables an opt-in initial Gram-Schmidt
Cholesky path for large matrices. The prototype performs the `ZPOTRF` step with
cuSOLVER and then computes `inv(U)` through cuBLAS `ZTRSM`, preserving the
existing CPU LAPACK fallback. The harness case is
`gpu_resident_stack_density_1cov_addoproj_cusolver_gram`; it keeps the
conservative default threshold at
`CPPAW_CUSOLVER_ACC_GRAM_CHOLESKY_MIN_N=4096`, because the 2048-band Si64 case
is compute-faster but still wall-time neutral after the extra matrix transfer.

Run directories on Spark GB10:

```
2048 repeats:
tests/profile/si64/runs/cusolver-gram-refresh-baseline-2048-20260601
tests/profile/si64/runs/cusolver-gram-refresh-gram-2048-20260601
tests/profile/si64/runs/cusolver-gram-refresh-gram-force-2048-20260601

4096 probe:
tests/profile/si64/runs/cusolver-gram-probe-baseline-4096-20260601
tests/profile/si64/runs/cusolver-gram-probe-gram-4096-20260601

Final harness case validation:
tests/profile/si64/runs/cusolver-gram-final-case-4096-20260601

Terok x86/A40 build and path smoke:
tests/profile/si64/runs/cusolver-gram-terok-smoke-2048-20260601
```

| Empty bands | Case | Repeats | Wall time | `PAW_GRAM_SOLVE` | `PAW_ETOT_SETUP_GRAM` | Transfer estimate | Energy check |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 2048 | ADDOPROJ stack | 3 | 30.27 s | 0.8527 s | 5.1160 s | 2.7307 GB | yes |
| 2048 | + cuSOLVER Gram | 3 | 30.22 s | 0.4436 s | 4.6776 s | 3.1853 GB | yes |
| 2048 | + cuSOLVER Gram + force DEDPRO | 3 | 30.31 s | 0.4491 s | 4.7266 s | 3.1416 GB | yes |
| 4096 | ADDOPROJ stack | 1 | 144.09 s | 5.9272 s | 17.6746 s | 7.0698 GB | yes |
| 4096 | + cuSOLVER Gram harness case | 1 | 136.61 s | 2.5734 s | 14.1987 s | 8.7827 GB | yes |

Representative 4096 rows show the intended split:
`CUSOLVER_ZPOTRF_GRAM=0.6185 s` and `CUBLAS_ZTRSM_GRAM=1.5685 s`, replacing
the CPU `LAPACK_ZPOTRF_GRAM=2.6692 s` and `LAPACK_ZTRTRI_GRAM=2.8986 s`.
The same code path also built on Terok x86/A40; a 2048-band smoke with the
threshold forced to `1` passed the Si64 energy check and recorded
`CUSOLVER_ZPOTRF_GRAM=0.1221 s`, `CUBLAS_ZTRSM_GRAM=0.2226 s`, and
`PAW_GRAM_SOLVE=0.4704 s`.
Conclusion: the path is worth keeping for large band/Gram cases, but should
remain opt-in and threshold-gated until a second large system or multi-step
case confirms the transfer trade-off.

## 2026-06-01 Ortho-X Resident DSYEVD And Gram Transform Copy Cleanup

The next residency cleanup keeps the real `WAVES_ORTHO_X` diagonalization on
the device when Ortho-X residency is active. Instead of calling the generic
copy-in/copy-out `LIB$DIAGR8` cuSOLVER route and then copying `U`/`EIG` back to
the GPU for the Newton iterations, `CUSOLVER_DSYEVD_PRESENT` symmetrizes
`CHIPSI` into resident `U`, runs cuSOLVER `DSYEVD`, and leaves `U` and `EIG`
resident for the following cuBLAS DGEMMs. The checked path is still available:
`CPPAW_CUSOLVER_ACC_CHECK=1` disables this present-device shortcut. The same
patch removes a redundant host `PSIINV=PSI` assignment in the Gram transform;
the existing explicit copy loop remains the single CPU/GPU copy path.

Run directories on Spark GB10:

```
512 smoke:
tests/profile/si64/runs/orthox-present-diag-smoke512-20260601

2048 repeats:
tests/profile/si64/runs/orthox-present-diag-2048-repeat-20260601

4096 probe:
tests/profile/si64/runs/orthox-present-diag-4096-20260601

Spark 4-rank smoke:
tests/profile/si64/runs/orthox-present-diag-4r-smoke512-20260601

Terok x86/A40 smoke:
tests/profile/si64/runs/orthox-present-diag-terok-smoke512-20260601
```

| Empty bands | Case | Repeats | Wall time | Transfer estimate | `PAW_ORTHO_X_DIAG` | `PAW_GRAM_TRANSFORM` | Energy check |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 4096 | previous cuSOLVER Gram harness case | 1 | 136.61 s | 8.7827 GB | 1.8607 s | 5.0584 s | yes |
| 4096 | + resident DSYEVD / no redundant `PSIINV` host copy | 1 | 137.85 s | 8.0690 GB | 1.7410 s | 4.8103 s | yes |
| 2048 | resident DSYEVD / no redundant `PSIINV` host copy | 3 | 30.89 s | 2.5413 GB | 0.2572 s | 1.2923 s | yes |

Spark also rebuilt both the serial and parallel residency-profile targets. The
4-rank smoke used the non-density `gpu_resident_addoproj` case, passed the
energy check, and recorded `CUSOLVER_DSYEVD_PRESENT` in rank-local profiles.
Terok rebuilt the serial target and passed the 512-band density-stack smoke
with `CUSOLVER_DSYEVD_PRESENT=0.0283 s`.

Conclusion: this is a memory-traffic/residency cleanup, not yet a wall-time
win. It removes the 4096-band `ACC_COPY_CUSOLVER_DSYEVD` traffic and trims the
Gram transform envelope, but the dominant cost remains the Ortho-X DGEMM
iteration sequence and force phase. Keep using it as enabling work for broader
wavefunction residency rather than as a standalone speedup claim.

## Recommended Next Benchmark

Use the focused default comparison for routine checks:

```
cd tests/profile/si64
TEST=si64_bands EMPTY_BANDS=1024 NSTEPS=1 RUN_GPU_ALL=no ./run_nvhpc_standard.sh
```

`NSTEPS=1` is now the harness default for this standard comparison, so the
explicit setting above is mostly documentation for reproduced Spark refreshes.

Use the full diagnostic sweep only when comparing library combinations:

```
cd tests/profile/si64
./run_gpu_library_matrix.sh
```

The library-matrix wrapper delegates to `run_nvhpc_standard.sh` with the same
1024-band, one-step shape and expands the GPU case list to the all-library,
single-library, library-disabled, NVLAMATH/NVBLAS, and memory-mode diagnostics.

Use the larger resource comparison when checking whether one MPI rank plus one
GPU still beats one-rank and eight-rank CPU/NVHPC references for the
band/orthogonalization-heavy path:

```
cd tests/profile/si64
./run_gpu_resource_comparison.sh
```

This wrapper codifies the `EMPTY_BANDS=2048`, `NSTEPS=1` focused case used in
the 2026-06-01 Spark/Terok comparison above.
