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

Set `CPPAW_CUFFT_ACC=1` to enable the native path. The runtime default is
conservative: only batches with at least `CPPAW_CUFFT_ACC_MIN_ELEMENTS=1000000`
elements are offloaded unless the environment overrides the threshold. Set the
threshold to `0` only for force-all diagnostics of the small-FFT overhead.

To profile the combined native GPU paths on one GPU, build an
`nvhpc_gpu_acc_*` target. This enables explicit cuBLAS by default, keeps native
cuFFT opt-in, and uses cuSOLVER only above its default size threshold. The same
binary can selectively force or disable each accelerator path:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_residency_profile
cd tests/profile/si64
NSTEPS=1 CASES="cpu nvhpc_cpu gpu_resident gpu_off" ./run_benchmark.sh
NSTEPS=1 RANKS=8 CASES="cpu nvhpc_cpu" ./run_benchmark.sh
```

`run_benchmark.sh` writes `benchmark.tsv` and `benchmark.md` for each run. The
summary includes both wall time and rank-normalized wall time (`wall_rank_s`),
the primary instrumented rank-seconds (`rank_s`), and a residual
`gap_s = wall_rank_s - rank_s`. Use `gap_s` and `coverage_pct` to decide
whether the current CSV timers already explain the run or whether additional
instrumentation is needed. Diagnostic Plane-wave FFT local/MPI-envelope timers
are reported separately as `pw_trace_s`; the GTOR and RTOG subsets are also
reported as `pw_gtor_s` and `pw_rtog_s`. These columns are intentionally kept
out of `rank_s` because they subdivide the existing `PW_FFT_*_TOTAL` envelope.
High-level `PHASE_*` timers are reported separately as `phase_s` and
`phase_gap_s = wall_rank_s - phase_s`; they are also kept out of `rank_s`
because they are coarse envelopes around existing numerical kernel timers. Use
the phase columns to localize unexplained wall time before adding lower-level
kernel instrumentation. The aggregate copy estimate `copy_gb` is split into
semantic buckets for the GPU-residency work: `copy_wave_gb` for wavefunction
arrays, `copy_proj_gb` for projector/projection arrays, `copy_offden_gb` for
off-site density-matrix transfers, and `copy_denmat_gb` for one-center
density-matrix transfers. `transfer_gb` is the combined host/device movement
estimate (`copy_gb + update_gb`). OpenACC `ACC_UPDATE_*` rows are reported
separately as `update_gb` with matching `update_wave_gb`, `update_proj_gb`,
`update_offden_gb`, and `update_denmat_gb` buckets, so required boundary
refreshes are visible without inflating the pure copy estimate. These buckets
are subsets of `copy_gb` or `update_gb` and are meant to show whether a
residency change actually removes the expected data motion. `profile_summary.py`
uses the same buckets in its category summary and prints their GB totals for
single-run CSV inspection.
For the bundled `si64` and `si64_bands` cases, the harness also checks the
final constant energy against the built-in reference (`EXPECTED_ENERGY`,
default `302.280854`) with `ENERGY_TOL=1e-5`. A run with normal termination but
an energy mismatch is reported as `ok=no` and carries `energy_delta` in the TSV
and Markdown summaries.
The harness also records inherited `CPPAW_GPU_*`, `CPPAW_CUBLAS_ACC_*`,
`CPPAW_CUSOLVER_ACC_*`, `CPPAW_CUFFT_ACC*`, and `CPPAW_GRAM_CHOLESKY`
environment switches in `run.env` and the benchmark `env` column. Case-specific
settings are appended after inherited settings, so explicit `CASES` keywords
remain reproducible and override broader shell defaults.

`WAVES$ETOT` is split further by `PAW_ETOT_*` rows, and the initial
Gram-Schmidt setup is split by `PAW_GRAM_*` rows. These are nested PAW
diagnostic envelopes: use them to identify the next target, not as additive
wall-clock accounting.
The real safe-orthogonalization solver is split further by `PAW_ORTHO_X_DIAG`,
`PAW_ORTHO_X_RESIDUAL`, `PAW_ORTHO_X_UPDATE`, and
`PAW_ORTHO_X_ITERATIONS`. `PAW_ORTHO_X_RESIDUAL_MATMUL` and
`PAW_ORTHO_X_RESIDUAL_CHECK` subdivide the residual row, while
`PAW_ORTHO_X_UPDATE_TRANSFORM`, `PAW_ORTHO_X_UPDATE_SCALE`,
`PAW_ORTHO_X_UPDATE_BACKTRANSFORM`, `PAW_ORTHO_X_UPDATE_APPLY`, and
`PAW_ORTHO_X_UPDATE_SYM` subdivide the update row. These rows are nested inside
`PAW_ORTHO_SOLVE`. In residency-profile builds, the broader `WAVES_ORTHO_X`
iteration-workspace residency is enabled by default and its
direct present-device cuBLAS calls are also reported as
`CUBLAS_DGEMM_ORTHOX_RESIDUAL`, `CUBLAS_DGEMM_ORTHOX_TRANSFORM`, and
`CUBLAS_DGEMM_ORTHOX_BACKTRANS`.
Residency-profile builds enable the initial Gram-Schmidt Cholesky solve by
default. It replaces only the initial `WAVES$GRAMMSCHMIDT` solve with a LAPACK
Cholesky orthogonalization when the overlap matrix is positive definite;
otherwise it falls back to the legacy solver. Set `CPPAW_GRAM_CHOLESKY=0` or
use `gpu_resident_gram_legacy` to compare against the old path; use
`gpu_resident_gram_cholesky` to force the new path explicitly.

The `nvhpc_gpu_acc_residency_*` target keeps the same accelerator choices but
defaults to `CPPAW_GPU_RESIDENCY=1`. Set `CPPAW_GPU_RESIDENCY=0` to disable the
resident mode in the same binary. This currently switches the
cuBLAS scalarproduct copy wrapper to `present_or_copyin`, so projection loops
can reuse wavefunction arrays already held by an outer OpenACC data region. For
non-superwave projections where all atom blocks pass the cuBLAS threshold, it
also scatters the projection result into `PROPSI` on the device and copies the
final projection array back once. The default residency profile now builds and
caches the full per-atom projector block `PRO` on the GPU from resident
`GSET%PRO`, `GSET%YLM`, and the updated structure factors, instead of expanding
`PRO` on the host and copying it for every atom/projection use.
`WAVES_PROJECTIONS` and eligible `WAVES_ADDPRO` calls share this cache while the
geometry, grid id and projector dimensions stay unchanged. Set
`CPPAW_GPU_PRO_EXPANSION=0` (or the longer alias
`CPPAW_CUBLAS_ACC_PRO_EXPANSION=0`) to compare against the host-expansion path;
the benchmark case is `gpu_resident_pro_host`. Set `CPPAW_GPU_ADDPRO_CACHE=0`
to keep the GPU projection cache but route `WAVES_ADDPRO` through the previous
host-expansion/addproduct path; the benchmark case is
`gpu_resident_addpro_host`. The narrower `CPPAW_GPU_ADDPRO_CACHE_HPSI=0` and
`CPPAW_GPU_ADDPRO_CACHE_OPSI=0` switches isolate the Hamiltonian and overlap
wavefunction `WAVES_ADDPRO` contexts without disabling the shared projection
cache. The orthogonalization overlap
section keeps `PSIM`/`OPSI` resident across the projection and pseudo-overlap
calls, and the same mode routes eligible `WAVES_OVERLAP`
scalarproducts through a present-input cuBLAS wrapper; for inversion-symmetric
superwave overlaps it also keeps the `<PSI_+|PSI_+>` part on the same
present-input path and batches the `<PSI_-|PSI_+>` inversion pass into one larger
scalarproduct when the cuBLAS overlap threshold allows it. `WAVES_GRAMSCHMIDT`
also has a narrow resident projection/overlap region for wavefunctions outside
the main orthogonalization loop. The force loop keeps `THIS%PSI0` resident
across per-atom `WAVES_DEDPRO` MATMUL calls; when HPSI residency is also
enabled, that same `PSI0` device copy is carried forward to the following
`WAVES$HPSI`/expectation/Hamiltonian boundary. Set
`CPPAW_GPU_FORCE_PSI_RESIDENCY=0` to compare against the previous per-atom copy
behavior. The one-center overlap contraction packs the flattened cuBLAS input
matrices on the GPU and is enabled by default in residency-profile builds; set
`CPPAW_GPU_1COVERLAP=0` to compare against the host contraction path. The real
safe-orthogonalization loop also has an opt-in diagnostic
that keeps the constant `CHICHI` and `U` matrices resident across repeated
`LIB$MATMULR8` calls while leaving the host-updated temporary outputs on the
previous copy-back path; set `CPPAW_GPU_ORTHO_CONST_RESIDENCY=1` to test it. It
is disabled by default because Spark Si64 smokes reduced copy volume but did not
improve wall time. A broader default path keeps the real `WAVES_ORTHO_X`
iteration workspace (`LAMBDA`, `GAMN`, `HAUX`, and scalar loop inputs) resident
and routes the large transform pairs through present-input cuBLAS calls; set
`CPPAW_GPU_ORTHO_X_RESIDENCY=0` or use `gpu_resident_orthox_off` to compare
against the previous copy-heavy path. Set
`CPPAW_CUBLAS_ACC_INVERSION_BATCH=0` to keep the older per-column inversion
scalarproduct path for comparison. It also lets
`ZGEMM_NN` addproduct calls reuse a present output matrix, which targets
`WAVES_ADDPRO`. An opt-in diagnostic can keep the orthogonalization `OPSI`
wavefunction resident through the projection, overlap, and `WAVES_ADDOPSI`
phase on non-stress paths where all atom blocks pass the addproduct threshold;
set `CPPAW_GPU_OPSI_RESIDENCY=1` or use `gpu_resident_opsi` to test it. For
non-superwave paths, OPSI can be resident from its build and mass scaling. For
inversion-symmetric superwave paths, OPSI is deliberately staged only after the
host build and host mass scaling because the fully resident build/scale variant
was energy-invalid in Si64. It is disabled by default until broader benchmarks
show that the reduced copy volume also improves wall time. Another opt-in
diagnostic keeps `HPSI` resident after the Hamiltonian-side `WAVES_ADDPRO`
update and reuses it for the immediate expectation and full-Hamiltonian overlap
calls. Eligible paths also keep `PSI0` present across the same local
`WAVES$ETOT` HPSI/expectation/Hamiltonian boundary and delete that ETOT-local
`PSI0` residency before leaving the energy evaluation; set
`CPPAW_GPU_HPSI_RESIDENCY=1` or use `gpu_resident_hpsi` to test it. The longer
alias is `CPPAW_CUBLAS_ACC_HPSI_RESIDENCY`. When HPSI residency is active and
the input `PSI` is already present, the final `WAVES_VPSI` kinetic/bucket
G-space update can finish on the device. The host-side FFT/RTOG boundary still
requires one `HPSI` host-to-device transfer, recorded as
`ACC_COPY_VPSI_HPSI_IN`; if `HPSI` is already present on the device, the same
required transfer is recorded as `ACC_UPDATE_VPSI_HPSI_IN` and the device copy
is refreshed explicitly before the kinetic/bucket update. The following
`WAVES_ADDPRO` consumer should then record `ACC_PRESENT_HPSI_ADDPRO` instead of
`ACC_COPY_HPSI_ADDPRO_IN`. A
separate opt-in
diagnostic propagates `PSIM` on the GPU and immediately copies the updated
wavefunction back before the following host-side projection work; set
`CPPAW_GPU_PSIM_PROPAGATE=1` or use `gpu_psim_propagate` /
`gpu_hpsi_psim_propagate`. The older `CPPAW_GPU_PSIM_RESIDENCY` alias is still
accepted for compatibility, but true cross-orthogonalization PSIM residency needs
broader projection/PRO residency first. These switches are disabled by default
because the extra propagation inputs and output copy can outweigh the kernel
offload. Set `CPPAW_GPU_PSIM_PHASE_RESIDENCY=1` or use
`gpu_resident_psim_phase` / `gpu_resident_hpsi_psim_phase` to keep the updated
`PSIM` array present from `WAVES$PROPAGATE` into the immediately following
orthogonalization block and copy it back only when orthogonalization finishes.
This cross-phase path is still opt-in because it relies on the current timestep
ordering and needs broader benchmarks before promotion. The propagation kernel
uses `present_or_copy` for `PSIM` and `present_or_copyin` for `PSI0`/`HPSI`, so
the broader resident region can reuse already-present wavefunction data without
changing the kernel body.
Generic
resident cuBLAS wrappers split their copy accounting
into `ACC_PRESENT_CUBLAS_*` and `ACC_COPY_CUBLAS_*` rows so already-resident
inputs are counted separately from real transfer estimates. The overlap region
uses `ACC_PRESENT_ORTHO_*` / `ACC_COPY_ORTHO_*` rows for its outer wavefunction
arrays and `ACC_COPY_CUBLAS_ZSPROD_OVL_RES` for the per-call output copy.
Projector-residency diagnostics include `ACC_BUILD_PRO_CACHE`,
`ACC_PRESENT_PRO_CACHE_REUSE`, `ACC_PRESENT_PROJ_PRO_CACHE`,
`ACC_PRESENT_ADDPRO_<ctx>_CACHE`, `CUBLAS_ZGEMM_ADDPRO_CACHE`, and the
disappearance or reduction of `ACC_COPY_PROJ_PRO_IN`. The resident overlap
cuBLAS kernels are timed separately as `CUBLAS_ZHERK_OVL_RES`,
`CUBLAS_ZGEMM_OVL_RES`, and `CUBLAS_ZGEMM_OVL_RES_INV` for the superwave
inversion contribution. This is the recommended NVHPC GPU performance path for
the larger Si64 band benchmarks. Keep
`gpu_resident_nosync` and `gpu_resident_orthox_nosync` as diagnostic candidates
only; Spark Nsight traces show that removing the explicit post-cuBLAS
synchronization mostly shifts waiting time into later stream synchronizations or
copy calls for this workload. The residency mode keeps the orthogonalization
`PSIM`/`OPSI` wavefunction pair resident from the projection/overlap phase
through `WAVES_ADDOPSI`. The `PSIM` region is recorded as split
`ACC_COPY_ORTHO_PSIM_IN` / `ACC_COPY_ORTHO_PSIM_OUT` rows when it is not
already present, and the ADDOPSI update itself should report
`ACC_PRESENT_ADDOPSI_PSIM`. For non-inversion wave sets, `OPSI` and `LAMBDA`
are copied into the addproduct region. With the inversion-batch GPU path enabled,
the inversion-symmetric path keeps `OPSI` present, creates the inverted
`OPSIINV` on the device, and reports it as `ACC_PRESENT_ADDOPSI_OPSIINV_TINV`.

The residency profile also records semantic OpenACC present checks for the PAW
wavefunction arrays that dominate this follow-up. `ACC_PRESENT_*` rows count
places where an array was already resident, while matching `ACC_COPY_*` rows add
the estimated bytes for a required host/device transfer. The tracked arrays are
`PSIM`/`OPSI` in the orthogonalization region and `WAVES_ADDOPSI`,
with `OPSI`/`LAMBDA` tracked for the non-inversion data region, `PSI` and
`PROPSI` in `WAVES_PROJECTIONS`, and context-specific projection `PSI` rows
such as `ACC_COPY_PROJ_SETUP0_PSI_IN`, `ACC_COPY_PROJ_GRAM_PSI0_PSI_IN`, and
`ACC_COPY_PROJ_ORTHO_PSIM_PSI_IN`. `WAVES_ADDPRO` has context-specific
`PSI`/`PROPSI` rows (`HPSI` and `OPSI`). Superwave projections use the same
resident cuBLAS projection path as ordinary wavefunctions, with the factor-two
and gamma-point correction applied on device before copying `PROPSI` back for
the communicator combine. The ADDPRO `PSI` rows are split into `*_PSI_IN` and
`*_PSI_OUT` estimates because the projector addition updates the wavefunction.
They are emitted for both the resident projector-cache path and the
host-expansion fallback path; in the fallback path, per-atom `PRO`/`PROPSI`
transfers remain in the generic cuBLAS `ZGEMM_NN` copy rows to avoid double
counting. The
wavefunction overlap `ZSPROD` rows are also tagged by `WAVES_OVERLAP` caller
context, for example `ACC_COPY_ZSP_ORTH_PSIM_P1_IN`,
`ACC_PRESENT_ZSP_GRAM_PSI0_P1`, and `ACC_COPY_ZSP_HAMILTON_P2_IN`.
The suffixes `P1`, `P2`, `OUT`, and `OVL` identify the first input, second
input, generic scalarproduct output, and resident-overlap output rows. The
Gram-Schmidt setup keeps `PSI` resident through the final wavefunction transform
and records that outer input/output region with context-specific rows such as
`ACC_COPY_GRAM_PSI0_PSI_IN`, `ACC_COPY_GRAM_PSI0_PSI_OUT`,
`ACC_COPY_GRAM_PSIM_PSI_IN`, and `ACC_COPY_GRAM_PSIM_PSI_OUT`. The
transform scratch `PSIINV` is created on the device from resident `PSI`,
recorded as `ACC_PRESENT_GRAM_PSIINV`, while the transform matrices are tracked
as `ACC_COPY_GRAM_X*`. The generic cuBLAS scalarproduct and `ZGEMM_NN` wrappers
also use these rows for their residency paths. HPSI residency records the
ETOT-local `PSI0` boundary as `ACC_COPY_HPSI_PSI0_IN` when it must create the
device copy, and the immediate consumers as `ACC_PRESENT_EXPECT_PSI0` and
`ACC_PRESENT_HAMILTON_PSI0` when they can reuse that resident buffer. If the
force-loop residency created the `PSI0` device copy earlier in the same
`WAVES$ETOT`, `WAVES$HPSI` records `ACC_PRESENT_HPSI_PSI0` instead of
`ACC_COPY_HPSI_PSI0_IN`. The `WAVES_VPSI` finish path records resident input
`PSI` as `ACC_PRESENT_VPSI_PSI`, the still-required host FFT/RTOG output
boundary as `ACC_COPY_VPSI_HPSI_IN`, or as `ACC_UPDATE_VPSI_HPSI_IN` when an
existing HPSI device allocation must be refreshed, and small kinetic/bucket
inputs as `ACC_COPY_VPSI_G2_IN` or `ACC_COPY_VPSI_BUCKET_IN`. When this
producer-side
device finish is active, the following ADDPRO step records
`ACC_PRESENT_HPSI_ADDPRO`. The
one-center overlap cuBLAS
path splits its estimated transfers into `ACC_COPY_1COV_PROJ_IN`,
`ACC_COPY_1COV_MAT_OUT`, and, for inversion-symmetric superwave cases,
`ACC_COPY_1COV_CMAT_OUT`; use those rows to decide whether a future optimization
should target packed projector inputs or overlap-matrix outputs.
One-center density-matrix profiling uses `PAW_DENMAT_*` rows to split the
previous `PAW_ETOT_DENMAT` envelope into occupation setup, site setup, inner
density/energy loops, accumulation, MPI combine, and spin conversion. Off-site
density-matrix setup is split into `PAW_OFFDEN_*` rows, and the off-site
accumulation envelope is split further into `PAW_OFFDEN_SUM_SETUP`,
`PAW_OFFDEN_SUM_ZERO`, `PAW_OFFDEN_SUM_LOCAL`, and
`PAW_OFFDEN_SUM_COMBINE`. These rows are CPU-side
instrumentation for deciding whether a later GPU kernel should target
`WAVES_DENMAT` itself, the projection copy/setup edges, or off-site bookkeeping.
An opt-in diagnostic, `CPPAW_GPU_DENMAT_ENERGY=1`, offloads the
time-inversion `WAVES_DENMAT` energy/Lambda contraction with OpenACC and records
`ACC_KERNEL_DENMAT_ENERGY_TINV` plus `ACC_COPY_DENMAT_ENERGY_TINV`; the harness
cases are `gpu_resident_denmat_energy` and
`gpu_resident_hpsi_denmat_energy`. The outer `WAVES$DENMAT` setup records
`PAW_DENMAT_LAGR_SETUP`, builds `LAGR=LAMBDA*OCC` once per k-point/spin, and,
when the diagnostic is active, keeps that block resident across the atom loop.
The LAGR device lifetime is visible through `ACC_COPY_DENMAT_LAGR_IN` and
`ACC_PRESENT_DENMAT_LAGR`; the remaining per-site output copy is reported as
`ACC_COPY_DENMAT_ENERGY_TINV`. It is disabled by default because the current
Si64 wall time is still neutral even though the DENMAT envelope and transfer
estimate shrink.
The scalar time-inversion off-site DENMAT diagnostic is enabled with
`CPPAW_GPU_OFFDEN_LOCAL=1`; `CPPAW_OFFDEN_BLAS=1` is kept as a shorter alias.
It rewrites the `TINV`/`NDIM=1` off-site local contraction as packed
`ZGEMM` calls and records `PAW_OFFDEN_BLAS_PACK`,
`ZGEMM_OFFDEN_TINV_NDIM1`, and `PAW_OFFDEN_BLAS_ACCUM`. This is still a
host-data BLAS prototype rather than a full resident GPU path, so it is
disabled by default and should be used to quantify whether keeping projector
buffers resident on the GPU would be worthwhile.
`CPPAW_GPU_OFFDEN_CUBLAS=1` or `CPPAW_CUBLAS_ACC_OFFDEN=1` additionally
tries cuBLAS for that packed `ZGEMM(N,T)` through
`CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP`; the forced harness cases set this threshold
to 1. Spark Si64 shows that the per-neighbor cuBLAS prototype is slower inside
`PAW_OFFDEN_SUM_LOCAL` than host BLAS because the matrices are tiny and copied
for every neighbor. Keep it as a diagnostic only; a useful GPU version should
batch neighbors and/or keep packed projector buffers resident.
`CPPAW_GPU_OFFDEN_CUBLAS_BATCH=1` or `CPPAW_CUBLAS_ACC_OFFDEN_BATCH=1`
switches the cuBLAS diagnostic from one tiny GEMM per neighbor to a stacked
formulation: neighbors with the same first atom and second-projector size are
packed into one wider `ZGEMM(N,T)` call. The chunk length is controlled by
`CPPAW_GPU_OFFDEN_BATCH_SIZE` and defaults to 64. It records
`PAW_OFFDEN_BATCH_PACK`, `CUBLAS_ZGEMM_OFFDEN_TINV_STACK`, and
`PAW_OFFDEN_BATCH_ACCUM`. This is still disabled by default because the current
prototype keeps packing host buffers; use it to quantify whether a resident
projector/off-site backend is worth implementing.
`CPPAW_GPU_OFFDEN_DEVICE_PACK=1` or
`CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_PACK=1` moves that stacked packing step into an
OpenACC region: `PROJ` and `OCC` are copied once for the off-site pass, A/B are
packed on the GPU, cuBLAS uses present data, and only the stacked `WORK` block
is copied back. It records `PAW_OFFDEN_DEVICE_PACK`,
`CUBLAS_ZGEMM_OFFDEN_TINV_DPACK`, `ACC_COPY_OFFDEN_DPACK_PROJ_IN`,
`ACC_COPY_OFFDEN_DPACK_META_IN`, and `ACC_COPY_OFFDEN_DPACK_WORK_OUT`. Spark
Si64 shows this is useful for the 1 MPI rank / 1 GPU comparison, but can be a
negative diagnostic when several MPI ranks share one GPU.
`CPPAW_GPU_OFFDEN_DEVICE_ACCUM=1` or
`CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_ACCUM=1` keeps the device-pack path active and
accumulates the complex `WORK` blocks into a flat real off-site matrix buffer
on the GPU. It then copies that flat buffer back once per k-point/spin pass
instead of copying each batch's complex `WORK`, recording
`PAW_OFFDEN_DEVICE_ACCUM`, `ACC_COPY_OFFDEN_DPACK_FLAT_OUT`, and
`PAW_OFFDEN_FLAT_ACCUM_SCATTER`. This is an opt-in diagnostic for reducing host
round-trips before a fully device-resident off-site accumulation path exists.
`CPPAW_GPU_PROJ_RESIDENCY=1` or `CPPAW_CUBLAS_ACC_PROJ_RESIDENCY=1` keeps
`THIS%PROJ` present after `WAVES$PROJECTIONS` and the `K`-communicator combine.
That lets downstream `PRESENT_OR_COPYIN` users reuse projections instead of
copying them again. The profile rows are `ACC_COPY_THIS_PROJ_IN`,
`ACC_PRESENT_THIS_PROJ`, and, for off-site device packing,
`ACC_PRESENT_OFFDEN_DPACK_PROJ`.
Inversion-symmetric Hermitian/symmetric scalarproducts no longer include the unused second
wavefunction array in the OpenACC data region. These rows are meant to guide the
next change: extend resident regions only where the profile shows repeated
copies of the same wavefunction data.

For an all-library diagnostic binary, build `nvhpc_gpu_all_*`. This links NVPL
fallbacks, cuFFTW, native cuFFT/OpenACC, cuBLAS/OpenACC, cuSOLVER/OpenACC and
NVLAMATH into one executable. NVBLAS stays separate because it interposes BLAS
calls at run time rather than being an explicit kernel path. The `gpu_all`
benchmark case enables cuFFT with the same conservative threshold, while
`gpu_force_all` and `gpu_resident_force_all` still use
`CPPAW_CUFFT_ACC_MIN_ELEMENTS=0` by default for overhead diagnostics:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_all_profile
cd tests/profile/si64
NSTEPS=1 CASES="gpu gpu_all gpu_all_off" ./run_benchmark.sh
```

The all-library binary is intentionally not part of the default larger-band
preset. It is useful to verify that all optional NVIDIA libraries link and run
together, but small and medium Si64 band cases are dominated by cuFFTW/NVLAMATH
overhead in this configuration.

The `cpu` case uses the plain GNU/OpenBLAS/FFTW build (`profile` or
`profile_parallel`) as a pre-HPC-SDK reference. For MPI runs it defaults to the
system or `cppaw-gccmpi` `mpirun`; for parallel CPU runs, the harness first
tries to derive the matching MPI launcher and `LD_LIBRARY_PATH` entry from the
selected executable's resolved `libmpi`. Set `CPU_MPIRUN=...` or
`CPU_MPI_LIBDIR=...` to override either value on unusual installations.
Benchmark runs pin common
CPU threading variables to one thread by default (`OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `BLIS_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS`, `NVPL_NUM_THREADS`) and record those settings plus
linked BLAS/LAPACK/FFT libraries in each run directory's `metadata.txt`.

The Si64 benchmark harness uses these `CASES` keywords:

| Keyword | Meaning |
| --- | --- |
| `cpu` | Plain GNU/OpenBLAS/FFTW reference build. |
| `nvhpc_cpu` | NVIDIA HPC SDK CPU build. The actual BLAS/LAPACK/FFT backend is recorded through `ldd` in `metadata.txt`; on x86 this can be OpenBLAS/FFTW fallback rather than NVPL. |
| `nvpl` | Legacy alias for `nvhpc_cpu`; kept for old scripts and logs. |
| `nvblas` | Existing BLAS calls through NVIDIA's NVBLAS interposition layer. |
| `nvlamath` | NVIDIA HPC SDK NVLAMATH LAPACK/cuSOLVER wrapper path. |
| `cufftw` | cuFFTW wrapper for CP-PAW's existing FFTW3 calls. |
| `cufft` / `cufft_force_all` / `cufft_off` | Native cuFFT enabled with the conservative size threshold, forced for all FFTs, or disabled in the same binary. |
| `cublas` / `cublas_nosync` / `cublas_invbatch_off` / `cublas_off` | Explicit cuBLAS/OpenACC path with the default threshold, with the post-call device synchronization disabled, with inversion scalarproduct batching disabled, or disabled. |
| `cublas_conservative` | Explicit cuBLAS/OpenACC with a higher diagnostic threshold. |
| `cublas_projection_conservative` / `cublas_overlap_conservative` / `cublas_addproduct_conservative` / `cublas_matmul_conservative` | Explicit cuBLAS/OpenACC with only one kernel category raised to the conservative threshold. |
| `cusolver` / `cusolver_standard` / `cusolver_generalized` / `cusolver_off` | Explicit cuSOLVER/OpenACC forced for all eigensolvers, only standard eigensolvers, only generalized eigensolvers, or disabled. |
| `cusolver_conservative` | cuSOLVER/OpenACC with the production default size threshold. |
| `cusolver_generalized_conservative` | cuSOLVER/OpenACC with only generalized eigensolvers enabled at the production threshold. |
| `gpu` / `gpu_nosync` / `gpu_invbatch_off` | Combined GPU profile, with optional diagnostics that disable the cuBLAS post-call synchronization or inversion scalarproduct batching. |
| `gpu_projection_conservative` / `gpu_overlap_conservative` / `gpu_addproduct_conservative` / `gpu_matmul_conservative` | Combined GPU diagnostics with only one cuBLAS kernel category raised to the conservative threshold. |
| `gpu_resident` / `gpu_resident_nosync` / `gpu_resident_invbatch_off` | Recommended combined GPU profile with `CPPAW_GPU_RESIDENCY=1`; currently keeps selected wavefunction loops in OpenACC data regions for cuBLAS scalarproduct/projection/addproduct reuse, with diagnostics for synchronization and inversion batching. |
| `gpu_resident_no_cusolver` | Residency diagnostic with cuSOLVER disabled in the same residency binary. |
| `gpu_resident_pro_host` | Residency diagnostic with GPU projector expansion disabled via `CPPAW_GPU_PRO_EXPANSION=0`. |
| `gpu_resident_proj` | Opt-in diagnostic that keeps `THIS%PROJ` resident after projection setup via `CPPAW_GPU_PROJ_RESIDENCY=1`. |
| `gpu_resident_addpro_host` | Residency diagnostic with the GPU projection cache kept enabled but its `WAVES_ADDPRO` reuse disabled via `CPPAW_GPU_ADDPRO_CACHE=0`. |
| `gpu_resident_forcepsi_host` | Residency diagnostic with force-loop `THIS%PSI0` residency disabled via `CPPAW_GPU_FORCE_PSI_RESIDENCY=0`. |
| `gpu_resident_1coverlap` / `gpu_resident_1coverlap_host` | Residency diagnostics that force or disable the one-center overlap cuBLAS path via `CPPAW_GPU_1COVERLAP`. |
| `gpu_resident_orthoconst` | Residency diagnostic with opt-in `WAVES_ORTHO_X` constant-input residency enabled via `CPPAW_GPU_ORTHO_CONST_RESIDENCY=1`. |
| `gpu_resident_orthox` | Explicit residency default with the real `WAVES_ORTHO_X` iteration workspace kept on the GPU via `CPPAW_GPU_ORTHO_X_RESIDENCY=1`. |
| `gpu_resident_orthox_off` | Residency diagnostic that disables the `WAVES_ORTHO_X` iteration workspace residency via `CPPAW_GPU_ORTHO_X_RESIDENCY=0`. |
| `gpu_resident_orthox_nosync` | Diagnostic that combines `gpu_resident_orthox` with `CPPAW_CUBLAS_ACC_SYNC=0`; use for profiling synchronization overhead, not as the default. |
| `gpu_resident_opsi` | Opt-in residency diagnostic that keeps orthogonalization `OPSI` on the GPU through projection/overlap/`WAVES_ADDOPSI` via `CPPAW_GPU_OPSI_RESIDENCY=1`; superwave cases use conservative host build/scale staging before device residency. |
| `gpu_resident_hpsi` | Opt-in residency diagnostic that keeps `HPSI` and ETOT-local `PSI0` on the GPU from Hamiltonian-side `WAVES_ADDPRO` through the immediate expectation/Hamiltonian overlaps via `CPPAW_GPU_HPSI_RESIDENCY=1`. |
| `gpu_resident_denmat_energy` | Opt-in diagnostic that sets `CPPAW_GPU_DENMAT_ENERGY=1` and forces the time-inversion one-center DENMAT energy/Lambda OpenACC prototype for comparison. |
| `gpu_resident_hpsi_denmat_energy` | Combined diagnostic with HPSI residency and the DENMAT energy/Lambda OpenACC prototype enabled together. |
| `gpu_resident_offden_blas` | Opt-in diagnostic that sets `CPPAW_GPU_OFFDEN_LOCAL=1` and rewrites scalar `TINV` off-site DENMAT local work as packed BLAS. |
| `gpu_resident_hpsi_offden_blas` | Combined HPSI residency plus scalar `TINV` off-site DENMAT BLAS diagnostic. |
| `gpu_resident_denmat_energy_offden_blas` | Combined DENMAT energy/Lambda diagnostic plus scalar `TINV` off-site DENMAT BLAS diagnostic. |
| `gpu_resident_hpsi_denmat_energy_offden_blas` | Combined HPSI residency, DENMAT energy/Lambda diagnostic, and scalar `TINV` off-site DENMAT BLAS diagnostic. |
| `gpu_resident_offden_cublas` | Diagnostic that forces `CPPAW_GPU_OFFDEN_CUBLAS=1` and `CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP=1` for the scalar off-site DENMAT BLAS prototype. |
| `gpu_resident_hpsi_offden_cublas` | Combined HPSI residency plus scalar off-site DENMAT cuBLAS diagnostic. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas` | Combined HPSI residency, DENMAT energy/Lambda diagnostic, and scalar off-site DENMAT cuBLAS diagnostic. |
| `gpu_resident_hpsi_offden_cublas_batch` | Combined HPSI residency plus stacked scalar off-site DENMAT cuBLAS diagnostic. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_batch` | Combined HPSI residency, DENMAT energy/Lambda diagnostic, and stacked scalar off-site DENMAT cuBLAS diagnostic. |
| `gpu_resident_hpsi_offden_cublas_devicepack` | Combined HPSI residency plus stacked scalar off-site DENMAT cuBLAS with OpenACC device packing. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack` | Combined HPSI residency, DENMAT energy/Lambda diagnostic, and stacked scalar off-site DENMAT cuBLAS with OpenACC device packing. |
| `gpu_resident_hpsi_offden_cublas_devicepack_accum` | Device-pack off-site DENMAT diagnostic with GPU-side real-matrix accumulation. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_accum` | Combined DENMAT energy/device-pack diagnostic with GPU-side real-matrix accumulation. |
| `gpu_resident_hpsi_offden_cublas_devicepack_proj` | Device-pack off-site DENMAT diagnostic with HPSI and persistent `THIS%PROJ` residency enabled. |
| `gpu_resident_hpsi_offden_cublas_devicepack_proj_accum` | Device-pack off-site DENMAT diagnostic with HPSI, persistent `THIS%PROJ`, and GPU-side real-matrix accumulation. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj` | Combined DENMAT energy/device-pack off-site diagnostic with persistent `THIS%PROJ` residency enabled. |
| `gpu_resident_hpsi_denmat_energy_offden_cublas_devicepack_proj_accum` | Combined DENMAT energy/device-pack diagnostic with persistent `THIS%PROJ` and GPU-side real-matrix accumulation. |
| `gpu_psim_propagate` | Opt-in diagnostic that propagates `PSIM` on the GPU and copies it back before orthogonalization via `CPPAW_GPU_PSIM_PROPAGATE=1`. |
| `gpu_hpsi_psim_propagate` | Combined diagnostic with both `CPPAW_GPU_HPSI_RESIDENCY=1` and `CPPAW_GPU_PSIM_PROPAGATE=1`. |
| `gpu_resident_psim_phase` | Opt-in diagnostic that also sets `CPPAW_GPU_PSIM_PHASE_RESIDENCY=1`, leaving propagated `PSIM` present until orthogonalization copies it back. |
| `gpu_resident_hpsi_psim_phase` | Combined HPSI plus cross-phase PSIM propagation residency diagnostic. |
| `gpu_resident_hpsi_opsi` | Combined residency diagnostic with both `CPPAW_GPU_HPSI_RESIDENCY=1` and `CPPAW_GPU_OPSI_RESIDENCY=1`. |
| `gpu_resident_hpsi_opsi_proj` | Combined HPSI/OPSI diagnostic with persistent `THIS%PROJ` residency enabled. |
| `gpu_resident_hpsi_opsi_offden_cublas_devicepack_accum` | Combined HPSI/OPSI diagnostic with off-site DENMAT device packing and GPU-side real-matrix accumulation. |
| `gpu_resident_stack` | Focused residency stack keyword via `CPPAW_GPU_RESIDENCY_STACK=1`; enables the validated HPSI, OPSI, PROJ, DENMAT energy, and off-site device-pack accumulation combination. |
| `gpu_resident_stack_cufft` / `gpu_resident_stack_cufft_force` | Focused residency stack plus native cuFFT enabled with the conservative threshold, or forced for all `LIB$FFTC8` calls. |
| `gpu_resident_hpsi_opsi_denmat_energy_offden_cublas_devicepack_proj_accum` | Full residency diagnostic that combines HPSI, OPSI, DENMAT energy, persistent `THIS%PROJ`, and off-site device-pack accumulation. |
| `gpu_resident_projection_conservative` / `gpu_resident_overlap_conservative` / `gpu_resident_addproduct_conservative` / `gpu_resident_matmul_conservative` | Residency diagnostics with only one cuBLAS kernel category raised to the conservative threshold. |
| `gpu_resident_force_all` | Residency diagnostic that also forces cuFFT and small cuSOLVER offload. |
| `gpu_resident_off` | Residency binary with native cuFFT/cuBLAS/cuSOLVER disabled for same-executable fallback comparison. |
| `gpu_all` / `gpu_all_nosync` / `gpu_all_invbatch_off` | All-library GPU diagnostic build with cuFFTW/NVLAMATH linked and native cuFFT/cuBLAS/cuSOLVER enabled at run time; cuFFT uses the conservative threshold by default, and inversion scalarproduct batching can be disabled for comparison. |
| `gpu_all_3dfft` | All-library diagnostic build with the opt-in native cuFFT 3-D wrapper enabled as well, also threshold-gated by default. |
| `gpu_all_off` | Same all-library binary with native cuFFT/cuBLAS/cuSOLVER disabled; cuFFTW/NVLAMATH remain compiled in. |
| `gpu_force_all` | Diagnostic combined profile that forces cuFFT, cuBLAS and small cuSOLVER offload. |
| `gpu_3dfft` | Diagnostic combined profile that also enables the opt-in native cuFFT 3-D wrapper. |
| `gpu_managed` / `gpu_unified` | Separate combined GPU binaries built with NVHPC `-gpu=mem:managed` or `-gpu=mem:unified` for memory-residency experiments. |
| `gpu_no_cufft` / `gpu_no_cublas` / `gpu_no_cusolver` | Diagnostic ablations from `gpu_force_all`. |
| `gpu_off` | Same combined binary with all native GPU paths disabled. |

By default, `run_benchmark.sh` uses `RANKS=1` and
`CASES="cpu nvhpc_cpu gpu_resident gpu_off"`. This matches the small Si64
recommendation: compare one MPI rank with one GPU residency path against
one-rank CPU references and the same combined binary with native GPU paths
disabled. The larger follow-up and band sweeps add `gpu_resident_stack` by
default so HPSI, OPSI, persistent projection, DENMAT-energy, and off-site
device-pack residency are measured in the CPU/GPU resource comparison. Use
explicit `CASES=...` or `GPU_CASES=...` for diagnostic sweeps.

For a larger band/orthogonalization smoke test, use the `si64_bands` control
file with the same Si64 structure and more empty bands. The harness copies
`si64.strc` automatically when a variant-specific structure file is not present:

```
cd tests/profile/si64
TEST=si64_bands EMPTY_BANDS=128 NSTEPS=1 ./run_benchmark.sh
```

The stacked follow-up benchmark compares the larger-band case across the
resource split we want for the next optimization pass: one rank on one GPU,
one-rank CPU references, and eight-rank CPU/NVHPC references. It defaults to
`NSTEPS_LIST="1 3 10"` and accepts `EMPTY_BANDS_LIST` for a band-size sweep. It
also writes `combined_transfer_rows.md` and `combined_present_rows.md` so the
remaining data motion and resident reuse sites can be compared across the sweep:

```
cd tests/profile/si64
./run_followup.sh
```

For the current standard NVHPC comparison used before opening follow-up PRs,
run:

```
cd tests/profile/si64
./run_nvhpc_standard.sh
```

It defaults to `TEST=si64_bands`, `EMPTY_BANDS=1024`, `NSTEPS=3` and by default
compares the focused `gpu_resident*` paths on one GPU rank, including the
opt-in Ortho-X workspace-residency diagnostic, plus one-rank CPU and eight-rank
CPU/NVHPC references. Override `GPU_CASES`, `CPU_CASES`,
`EMPTY_BANDS`, `NSTEPS`,
`GPU_RANKS` or `CPU_RANKS` for a targeted sweep. It writes both
`combined_benchmark.md` and `combined_compare.md` so PR comments can include
raw timings and the one-GPU versus CPU-resource speedup view. It also writes
`combined_transfer_rows.md`/`.tsv` for the largest `ACC_COPY` and `ACC_UPDATE`
rows and `combined_present_rows.md`/`.tsv` for the most frequent
`ACC_PRESENT` rows across the suites, so the next residency change can be picked
from the measured remaining data-motion and reuse sites. Set `PROFILE_ROW_TOP`
or `PRESENT_ROW_TOP` to widen those reports.

Set `RUN_GPU_ALL=yes` to include the all-library cases `gpu_all` and
`gpu_all_off`. The default is `RUN_GPU_ALL=no` because the Spark Si64 matrix
showed the all-library path is useful as a diagnostic, not as a recommended
default.
Set `AUTO_BUILD_TARGETS=yes` (with `AUTO_BUILD_JOBS`) to automatically build all
required profile binaries before benchmarking.

The Spark C86C Si64 decision table is kept in
`tests/profile/si64/nvhpc_spark_benchmark_summary.md`.

For the focused off-site DENMAT residency comparison, run:

```
cd tests/profile/si64
./run_offden_focus.sh
```

It compares device-pack, flat device-accumulation, projection residency, and the
combined projection-plus-accumulation paths for both HPSI-only and
DENMAT-energy cases. By default it runs `EMPTY_BANDS=2048` on one GPU rank and
`SHARED_EMPTY_BANDS=512` on four ranks sharing the GPU. Override
`OFFDEN_CASES`, `GPU_RANKS`, `SHARED_GPU_RANKS`, `EMPTY_BANDS`, or
`RUN_SHARED_GPU=no` for narrower checks. It writes combined benchmark and
comparison Markdown so the device-pack variants can be compared against the
first successful focus-case baseline.

For the focused PSIM propagation comparison, run:

```
cd tests/profile/si64
./run_psim_focus.sh
```

It compares the default residency path, PSIM propagation, cross-phase PSIM
propagation residency, HPSI residency, and the combined HPSI-plus-PSIM
diagnostics. By default it runs
`EMPTY_BANDS=512` with one GPU rank and `SHARED_EMPTY_BANDS=512` with four
ranks sharing the GPU, then writes a combined TSV/Markdown summary next to the
run directories, including a comparison table against the focus baseline. Set
`RUN_LARGE_GPU=yes` to add a one-rank
`LARGE_EMPTY_BANDS=2048` sweep, or override `PSIM_CASES`, `GPU_RANKS`,
`SHARED_GPU_RANKS`, `EMPTY_BANDS`, and `RUN_SHARED_GPU=no` for narrower
checks.

For the PSIM lifecycle comparison across at least two time steps, run:

```
cd tests/profile/si64
./run_psim_lifecycle.sh
```

It uses `NSTEPS=2` by default and compares the same PSIM/HPSI diagnostics as
`run_psim_focus.sh`, but keeps the run shape centered on one GPU rank. This is
intended to expose copies around the propagation, orthogonalization, and next
time-step boundaries before broader cross-step wavefunction residency is enabled.
Because the one-step Si64 reference energy is not valid for multi-step dynamics,
the fixed energy check is disabled by default for this harness; compare the
reported energies between cases instead. It writes the normal combined benchmark
TSV/Markdown plus per-case `ACC_COPY` row summaries from `profile_copy_rows.py`;
the comparison Markdown reports each GPU-only focus case relative to the
selected group baseline. The helper also accepts `--op-prefix`, `--op-regex`,
and `--include-zero` for broader profile-row reports such as
present/update/timing rows; its size column is `gbyte`, so mixed copy, update,
and present reports do not imply every selected row is a copy.
Override `NSTEPS_LIST`, `EMPTY_BANDS_LIST`, `RUN_SHARED_GPU=yes`,
`RUN_LARGE_GPU=yes`, or `PSIM_LIFECYCLE_CASES` for wider sweeps.

For the VPSI/HPSI producer-boundary comparison, run:

```
cd tests/profile/si64
./run_vpsi_boundary.sh
```

It focuses on the current remaining HPSI boundary: `WAVES_VPSI` still receives
host-side FFT/RTOG output, then refreshes or creates the device-present `HPSI`
buffer for the following `WAVES_ADDPRO` consumer. The harness compares
`gpu_resident_hpsi`, `gpu_resident_hpsi_opsi`, `gpu_resident_stack`, and
stack-plus-cuFFT variants with `NSTEPS=2` by default and writes combined
benchmark tables plus selected profile tables for
`ACC_COPY`/`ACC_UPDATE`/`ACC_PRESENT` rows, separate seconds-sorted
`PAW_VPSI_*` timing rows, and separate seconds-sorted `PW_GTOR_*`/`PW_RTOG_*`
phase rows. Its comparison Markdown reports the stack and cuFFT variants
relative to the selected focus baseline.
Because multi-step Si64 energies are not the one-step reference, fixed energy
checking is disabled by default; compare energies between cases. Override
`VPSI_BOUNDARY_CASES`, `NSTEPS_LIST`, `EMPTY_BANDS_LIST`,
`RUN_SHARED_GPU=yes`, or `RUN_LARGE_GPU=yes` for wider sweeps.

For the larger orthogonalization preset used in the residency follow-up, run:

```
cd tests/profile/si64
./run_large_bands.sh
```

It defaults to `EMPTY_BANDS_LIST="128 256 512 1024"`, `NSTEPS_LIST=1` and
compares one-rank GPU residency, one-rank CPU, and eight-rank CPU/NVHPC
references. It inherits the follow-up benchmark's combined comparison,
transfer-row, and present-row reports. Set `RUN_GPU_ALL=yes` to include the
all-library diagnostic cases `gpu_all` and `gpu_all_off`.

For the longer validation preset used before promoting a residency diagnostic
to a default, run:

```
cd tests/profile/si64
./run_large_bands_long.sh
```

It defaults to `EMPTY_BANDS_LIST="512 1024"` and `NSTEPS_LIST="3 10"` and
compares `gpu_resident`, `gpu_resident_nosync`, `gpu_off`, one-rank CPU and
eight-rank CPU/NVHPC references. It inherits the follow-up benchmark's combined
comparison, transfer-row, and present-row reports.

For a short broad GPU exploration suite that includes the opt-in 3-D cuFFT path
and captures available NVIDIA libraries plus CUDA-aware MPI hints:

```
cd tests/profile/si64
NSTEPS=1 ./run_gpu_exploration.sh
```

The exploration run writes the same combined benchmark, comparison, transfer-row,
and present-row reports as the standard benchmark, plus `gpu_capabilities.txt`
when the capability helper is available.

To compare NVHPC memory modes, build the optional profile binaries and add the
cases explicitly:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_managed_profile
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_gpu_acc_unified_profile
cd tests/profile/si64
NSTEPS=1 CASES="gpu gpu_managed gpu_unified gpu_off" ./run_benchmark.sh
```

To test the explicit cuBLAS/OpenACC path for large complex `ZGEMM`/`ZHERK`
kernels, build an `nvhpc_cublas_acc_*` target. The default offload threshold is
`CPPAW_CUBLAS_ACC_MINFLOP=1e7`, which includes the projection GEMMs and was the
best Si64 threshold in the Spark C86C night run. Set `CPPAW_CUBLAS_ACC=0` to run
the same binary with the CPU fallback. Set `CPPAW_CUBLAS_ACC_SYNC=0` only
for diagnostic runs that compare the cost of the explicit device synchronization:

The global threshold is still the default for all cuBLAS call sites, but it can
be overridden by kernel category:

- `CPPAW_CUBLAS_ACC_PROJECTION_MINFLOP`: projector GEMMs in
  `WAVES_PROJECTIONS`
- `CPPAW_CUBLAS_ACC_OVERLAP_MINFLOP`: wavefunction overlap/orthogonalization
  scalar products
- `CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP`: additive projector/product updates
- `CPPAW_CUBLAS_ACC_MATMUL_MINFLOP`: generic library `MATMUL` replacements
- `CPPAW_CUBLAS_ACC_INVERSION_BATCH`: keep enabled by default to turn
  inversion-symmetry scalarproducts from many per-column cuBLAS calls into one
  batched scalarproduct; set to `0` for the previous path.
- `CPPAW_GPU_RESIDENCY_STACK`: disabled by default. Set to `1` to enable the
  focused Spark-validated residency stack in one switch: base residency,
  projector expansion/cache, HPSI and OPSI residency, persistent `THIS%PROJ`,
  one-center overlap, DENMAT energy offload, and off-site DENMAT cuBLAS
  device-pack accumulation. It also lowers the projection, overlap, addproduct,
  DENMAT, and off-site DENMAT thresholds for that run. Specific `CPPAW_GPU_*`
  or `CPPAW_CUBLAS_ACC_*` switches still override the corresponding part of
  the stack; disabling an off-site parent switch also disables its dependent
  child paths unless a later, more specific child switch explicitly re-enables
  the required parent path. The
  compatibility alias is `CPPAW_CUBLAS_ACC_RESIDENCY_STACK`.
- `CPPAW_GPU_PRO_EXPANSION`: keep enabled by default in residency-profile builds
  so GPU-resident `PRO` blocks are built once and reused by `WAVES_PROJECTIONS`
  and eligible `WAVES_ADDPRO` calls; set to `0` for the previous host-expansion
  path.
- `CPPAW_GPU_ADDPRO_CACHE`: keep enabled by default in residency-profile builds
  so `WAVES_ADDPRO` reuses the GPU-resident `PRO` cache; set to `0` to test
  projection caching without the cached addproduct path.
- `CPPAW_GPU_ADDPRO_CACHE_HPSI` and `CPPAW_GPU_ADDPRO_CACHE_OPSI`: keep enabled
  by default and override only the `HPSI` or `OPSI` `WAVES_ADDPRO` context. The
  harness cases are `gpu_resident_addpro_hpsi_host`,
  `gpu_resident_addpro_opsi_host`, and `gpu_resident_opsi_addpro_host`.
- `CPPAW_GPU_FORCE_PSI_RESIDENCY`: keep enabled by default in residency-profile
  builds so `WAVES$FORCE` reuses `THIS%PSI0` across the per-atom
  `WAVES_DEDPRO` MATMUL calls. When `CPPAW_GPU_HPSI_RESIDENCY=1` is also set,
  the same `PSI0` device copy is kept for the immediately following HPSI
  overlap consumers; set to `0` for the previous per-call copy path.
- `CPPAW_GPU_1COVERLAP`: keep enabled by default in residency-profile builds so
  `WAVES_1COVERLAP` uses the GPU-pack/cuBLAS contraction path; set to `0` for
  the host contraction path.
- `CPPAW_GPU_ORTHO_CONST_RESIDENCY`: disabled by default. Set to `1` to let
  `WAVES_ORTHO_X` reuse constant `CHICHI` and `U` inputs across repeated real
  MATMUL calls while keeping temporary outputs on the host-synchronized path.
- `CPPAW_GPU_ORTHO_X_RESIDENCY`: enabled by default in residency-profile builds
  so the real `WAVES_ORTHO_X` iteration workspace stays on the GPU and the large
  residual/transform/backtransform pairs use present-input cuBLAS calls; set to
  `0` for the previous copy-heavy path.
- `CPPAW_GPU_OPSI_RESIDENCY`: disabled by default. Set to `1` to keep
  orthogonalization `OPSI` resident through projection, overlap, and
  `WAVES_ADDOPSI` on eligible non-stress paths. Eligible paths keep OPSI
  resident from build and mass scaling. Inversion-symmetric superwave paths are
  only admitted when the projection threshold also selects the resident cuBLAS
  projection path, so the post-mass host snapshot is avoided and the later
  projection, overlap, and `WAVES_ADDOPSI` consumers use present device data.
  The
  compatibility alias is `CPPAW_CUBLAS_ACC_OPSI_RESIDENCY`.
- `CPPAW_GPU_HPSI_RESIDENCY`: disabled by default. Set to `1` to keep the
  Hamiltonian-side `HPSI` and matching `PSI0` inputs resident through the
  immediate expectation and full-Hamiltonian overlap calls within
  `WAVES$ETOT`. The compatibility alias is
  `CPPAW_CUBLAS_ACC_HPSI_RESIDENCY`.
- `CPPAW_GPU_PSIM_PROPAGATE`: disabled by default. Set to `1` to run
  `WAVES$PROPAGATE` on the GPU for non-stress steps and copy the updated `PSIM`
  back before orthogonalization. The compatibility aliases are
  `CPPAW_GPU_PSIM_RESIDENCY` and `CPPAW_CUBLAS_ACC_PSIM_RESIDENCY`; they do not
  imply true cross-orthogonalization residency in the current implementation.
- `CPPAW_GPU_PSIM_PHASE_RESIDENCY`: disabled by default and requires
  `CPPAW_GPU_PSIM_PROPAGATE=1`. Set to `1` to leave the propagated `PSIM` array
  present until the immediately following orthogonalization block copies the
  final orthogonalized wavefunction back. The compatibility aliases are
  `CPPAW_GPU_PSIM_KEEP_RESIDENT` and
  `CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY`.
- `CPPAW_GPU_DENMAT_ENERGY`: disabled by default. Set to `1` to offload the
  time-inversion one-center DENMAT energy/Lambda contraction with OpenACC. The
  compatibility alias is `CPPAW_CUBLAS_ACC_DENMAT_ENERGY`; use
  `CPPAW_GPU_DENMAT_MINFLOP` or `CPPAW_CUBLAS_ACC_DENMAT_MINFLOP` to adjust the
  offload threshold. Residency-profile builds precompute the corresponding
  `LAGR=LAMBDA*OCC` block once per k-point/spin and keep it on the GPU while
  the DENMAT diagnostic is active.

The benchmark harness exposes conservative diagnostic cases such as
`gpu_resident_projection_conservative`, `gpu_resident_overlap_conservative`,
`gpu_resident_addproduct_conservative`, and `gpu_resident_matmul_conservative`.
They keep the recommended `1e7` default for the other categories and raise only
one category to `CPPAW_CUBLAS_CONSERVATIVE_MINFLOP` (default `1e8`), unless a
category-specific conservative value is set.

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
case, or set `CPPAW_CUSOLVER_ACC=0` to run the same binary with the CPU/NVHPC
fallback. The global threshold can be split into
`CPPAW_CUSOLVER_ACC_STANDARD_MIN_N` for `DSYEVD/ZHEEVD` and
`CPPAW_CUSOLVER_ACC_GENERALIZED_MIN_N` for `DSYGVD/ZHEGVD`; exact overrides
`CPPAW_CUSOLVER_ACC_DSYEVD_MIN_N`, `CPPAW_CUSOLVER_ACC_ZHEEVD_MIN_N`,
`CPPAW_CUSOLVER_ACC_DSYGVD_MIN_N` and `CPPAW_CUSOLVER_ACC_ZHEGVD_MIN_N` are
also accepted. Use `cusolver_generalized` to force only the generalized path
and `cusolver_generalized_conservative` for the threshold-gated variant:

```
CPPAW_TOOLCHAIN=nvhpc src/Buildtools/paw_build.sh -c nvhpc_cusolver_acc_profile
cd tests/profile/si64
NSTEPS=1 RANKS=1 CASES="nvhpc_cpu cusolver cusolver_generalized cusolver_off" ./run_benchmark.sh
```

For a more targeted cuSOLVER/LAPACK follow-up on larger band matrices:

```
cd tests/profile/si64
./run_cusolver_focus.sh
```

It defaults to `EMPTY_BANDS_LIST="128 256 512"` and compares
`cusolver`, `cusolver_generalized`, `cusolver_generalized_conservative`,
`cusolver_off`, one-rank CPU/NVHPC and eight-rank CPU/NVHPC references. It
writes `combined_benchmark.md`, `combined_compare.md`,
`combined_solver_rows.md`, `combined_transfer_rows.md`, and
`combined_present_rows.md`. The solver-row report is sorted by seconds over
`LAPACK*` and `CUSOLVER*` rows, so the focus run shows directly which
eigensolver path was active and how much solver time remains.

For reproducible comparisons, use the benchmark harness:

```
cd tests/profile/si64
NSTEPS=20 RANKS=4 REPEATS=3 CASES="nvhpc_cpu cublas cublas_off" ./run_benchmark.sh
```

The harness creates timestamped directories under `tests/profile/si64/runs`,
writes per-run logs and profile CSV files, and emits a `benchmark.tsv` summary
with wall time, instrumented rank-seconds, category timings, final energy and
the optional energy delta. It also writes a Markdown table (`benchmark.md` or
`combined_benchmark.md`) that can be pasted directly into pull request comments.
The combined GPU exploration and follow-up harnesses also write
`combined_compare.md`, which reports wall-time speedups against the matching
one-rank and eight-rank CPU baselines plus the GPU residency copy buckets.

For a short Nsight Systems trace of the recommended combined GPU profile binary:

```
cd tests/profile/si64
CASE=gpu_resident NSTEPS=1 ./run_nsys.sh
```

`run_nsys.sh` accepts the same GPU-oriented case names as the benchmark harness,
for example `gpu_resident`, `gpu_resident_orthox`, `gpu_all`, `gpu_no_cufft`,
`cublas`, `cusolver`, and `cufft`. It writes `nsys_case.env` into the run
directory so a trace can be matched to the executable and runtime switches. If
Nsight leaves a SQLite export next to the report, the harness also writes
`nsys_sql_summary.txt` with the top CUDA kernels, runtime calls, memcpy totals,
and synchronization totals.
For `RANKS>1` the Nsight harness wraps `mpirun` and writes one combined
`nsys_mpi.nsys-rep` report, which is robust for CP-PAW's current `MPI_ABORT(0)`
shutdown path. For runs that finalize MPI normally, `NSYS_MPI_MODE=per_rank`
writes one report per rank.

For an overnight comparison that combines a longer 4-rank run, rank scaling,
a cuBLAS offload-threshold sweep and a short Nsight trace:

```
cd tests/profile/si64
./run_overnight.sh
```

The top-level run directory is written to `runs/latest_overnight`; the combined
benchmark table is `combined_benchmark.tsv`. The overnight run also writes
`combined_benchmark.md`, `combined_compare.md`, `combined_transfer_rows.md` and
`combined_present_rows.md` so the long run ends with the same speedup,
data-motion, and resident-reuse summaries as the shorter benchmark harnesses.

Set `RUN_NVLAMATH=yes` to add a short NVLAMATH comparison suite,
`RUN_CUFFTW=yes` to add a short cuFFTW comparison suite, `RUN_CUFFT=yes` to add
a short native cuFFT comparison, `RUN_GPU_ACC=yes` to compare one-rank GPU
against eight-rank CPU/OpenBLAS and CPU/NVHPC references, or
`RUN_CUSOLVER=yes` to add the same one-rank cuSOLVER versus eight-rank CPU/NVHPC
resource comparison.

The overnight defaults use the recommended production-style cases:
`MAIN_CASES="nvhpc_cpu cublas cublas_off"`,
`SCALING_CASES="nvhpc_cpu cublas"`,
`GPU_ACC_CASES="cpu nvhpc_cpu gpu_resident gpu_resident_stack gpu_resident_nosync gpu gpu_off"`
and `THRESHOLDS="1e7"`. Set
`RUN_GPU_DIAGNOSTICS=yes` to add `gpu_nosync`, `gpu_resident_stack_cufft`,
`gpu_force_all` and the `gpu_no_*` ablation cases, set
`RUN_BAND_BENCHMARK=yes` to add the
`si64_bands` larger-band matrix over one-rank GPU, one-rank CPU and eight-rank
CPU references, set `BAND_EMPTY_BANDS_LIST="128 256 512"` for a size sweep, or
override any of these variables for a wider run.

The capability helper can be run standalone:

```
src/Tools/Scripts/paw_gpu_capabilities.sh
```

It reports CUDA devices, NVIDIA HPC SDK library presence for cuBLASLt,
cuSPARSE, cuTENSOR, cuDSS, NCCL and NVSHMEM, and a best-effort CUDA-aware MPI
hint. Those libraries are profiled as future candidates; they are not linked
into CP-PAW unless a concrete code path uses them.

For an active CUDA-aware MPI smoke test, use:

```
src/Tools/Scripts/paw_cuda_aware_mpi_probe.sh
```

It compiles a tiny MPI/OpenACC allreduce probe and passes device pointers to
MPI. A failure here means GPU-resident MPI communication should stay disabled
for production runs until the MPI stack is configured appropriately.
