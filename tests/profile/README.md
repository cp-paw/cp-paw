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
are reported separately as `pw_trace_s`; they are intentionally kept out of
`rank_s` because they subdivide the existing `PW_FFT_*_TOTAL` envelope.
High-level `PHASE_*` timers are reported separately as `phase_s` and
`phase_gap_s = wall_rank_s - phase_s`; they are also kept out of `rank_s`
because they are coarse envelopes around existing numerical kernel timers. Use
the phase columns to localize unexplained wall time before adding lower-level
kernel instrumentation.
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
across per-atom `WAVES_DEDPRO` MATMUL calls; set
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
calls; set `CPPAW_GPU_HPSI_RESIDENCY=1` or use `gpu_resident_hpsi` to test it.
The longer alias is `CPPAW_CUBLAS_ACC_HPSI_RESIDENCY`. A separate opt-in
diagnostic propagates `PSIM` on the GPU and immediately copies the updated
wavefunction back before the following host-side projection work; set
`CPPAW_GPU_PSIM_PROPAGATE=1` or use `gpu_psim_propagate` /
`gpu_hpsi_psim_propagate`. The older `CPPAW_GPU_PSIM_RESIDENCY` alias is still
accepted for compatibility, but true cross-orthogonalization PSIM residency needs
broader projection/PRO residency first. These switches are disabled by default
because the extra propagation inputs and output copy can outweigh the kernel
offload.
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
`PSI`/`PROPSI` rows (`HPSI` and `OPSI`). The ADDPRO `PSI` rows are split into
`*_PSI_IN` and `*_PSI_OUT` estimates because the projector addition updates the
wavefunction. They are emitted for both the resident projector-cache path and
the host-expansion fallback path; in the fallback path, per-atom `PRO`/`PROPSI`
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
also use these rows for their residency paths. The one-center overlap cuBLAS
path splits its estimated transfers into `ACC_COPY_1COV_PROJ_IN`,
`ACC_COPY_1COV_MAT_OUT`, and, for inversion-symmetric superwave cases,
`ACC_COPY_1COV_CMAT_OUT`; use those rows to decide whether a future optimization
should target packed projector inputs or overlap-matrix outputs.
One-center density-matrix profiling uses `PAW_DENMAT_*` rows to split the
previous `PAW_ETOT_DENMAT` envelope into occupation setup, site setup, inner
density/energy loops, accumulation, MPI combine, and spin conversion. Off-site
density-matrix setup is split into `PAW_OFFDEN_*` rows. These rows are CPU-side
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
| `gpu_resident_addpro_host` | Residency diagnostic with the GPU projection cache kept enabled but its `WAVES_ADDPRO` reuse disabled via `CPPAW_GPU_ADDPRO_CACHE=0`. |
| `gpu_resident_forcepsi_host` | Residency diagnostic with force-loop `THIS%PSI0` residency disabled via `CPPAW_GPU_FORCE_PSI_RESIDENCY=0`. |
| `gpu_resident_1coverlap` / `gpu_resident_1coverlap_host` | Residency diagnostics that force or disable the one-center overlap cuBLAS path via `CPPAW_GPU_1COVERLAP`. |
| `gpu_resident_orthoconst` | Residency diagnostic with opt-in `WAVES_ORTHO_X` constant-input residency enabled via `CPPAW_GPU_ORTHO_CONST_RESIDENCY=1`. |
| `gpu_resident_orthox` | Explicit residency default with the real `WAVES_ORTHO_X` iteration workspace kept on the GPU via `CPPAW_GPU_ORTHO_X_RESIDENCY=1`. |
| `gpu_resident_orthox_off` | Residency diagnostic that disables the `WAVES_ORTHO_X` iteration workspace residency via `CPPAW_GPU_ORTHO_X_RESIDENCY=0`. |
| `gpu_resident_orthox_nosync` | Diagnostic that combines `gpu_resident_orthox` with `CPPAW_CUBLAS_ACC_SYNC=0`; use for profiling synchronization overhead, not as the default. |
| `gpu_resident_opsi` | Opt-in residency diagnostic that keeps orthogonalization `OPSI` on the GPU through projection/overlap/`WAVES_ADDOPSI` via `CPPAW_GPU_OPSI_RESIDENCY=1`; superwave cases use conservative host build/scale staging before device residency. |
| `gpu_resident_hpsi` | Opt-in residency diagnostic that keeps `HPSI` on the GPU from Hamiltonian-side `WAVES_ADDPRO` through the immediate expectation/Hamiltonian overlaps via `CPPAW_GPU_HPSI_RESIDENCY=1`. |
| `gpu_resident_denmat_energy` | Opt-in diagnostic that sets `CPPAW_GPU_DENMAT_ENERGY=1` and forces the time-inversion one-center DENMAT energy/Lambda OpenACC prototype for comparison. |
| `gpu_resident_hpsi_denmat_energy` | Combined diagnostic with HPSI residency and the DENMAT energy/Lambda OpenACC prototype enabled together. |
| `gpu_psim_propagate` | Opt-in diagnostic that propagates `PSIM` on the GPU and copies it back before orthogonalization via `CPPAW_GPU_PSIM_PROPAGATE=1`. |
| `gpu_hpsi_psim_propagate` | Combined diagnostic with both `CPPAW_GPU_HPSI_RESIDENCY=1` and `CPPAW_GPU_PSIM_PROPAGATE=1`. |
| `gpu_resident_hpsi_opsi` | Combined residency diagnostic with both `CPPAW_GPU_HPSI_RESIDENCY=1` and `CPPAW_GPU_OPSI_RESIDENCY=1`. |
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
`CASES="cpu nvhpc_cpu gpu_resident gpu_off"`. This matches the current Si64
recommendation: compare one MPI rank with one GPU residency path against
one-rank CPU references and the same combined binary with native GPU paths
disabled. Use explicit `CASES=...` for diagnostic sweeps.

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
`NSTEPS_LIST="1 3 10"` and accepts `EMPTY_BANDS_LIST` for a band-size sweep:

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
`GPU_RANKS` or `CPU_RANKS` for a targeted sweep.

Set `RUN_GPU_ALL=yes` to include the all-library cases `gpu_all` and
`gpu_all_off`. The default is `RUN_GPU_ALL=no` because the Spark Si64 matrix
showed the all-library path is useful as a diagnostic, not as a recommended
default.
Set `AUTO_BUILD_TARGETS=yes` (with `AUTO_BUILD_JOBS`) to automatically build all
required profile binaries before benchmarking.

The Spark C86C Si64 decision table is kept in
`tests/profile/si64/nvhpc_spark_benchmark_summary.md`.

For the larger orthogonalization preset used in the residency follow-up, run:

```
cd tests/profile/si64
./run_large_bands.sh
```

It defaults to `EMPTY_BANDS_LIST="128 256 512 1024"`, `NSTEPS_LIST=1` and
compares one-rank GPU residency, one-rank CPU, and eight-rank CPU/NVHPC
references. Set `RUN_GPU_ALL=yes` to include the all-library diagnostic cases
`gpu_all` and `gpu_all_off`.

For the longer validation preset used before promoting a residency diagnostic
to a default, run:

```
cd tests/profile/si64
./run_large_bands_long.sh
```

It defaults to `EMPTY_BANDS_LIST="512 1024"` and `NSTEPS_LIST="3 10"` and
compares `gpu_resident`, `gpu_resident_nosync`, `gpu_off`, one-rank CPU and
eight-rank CPU/NVHPC references.

For a short broad GPU exploration suite that includes the opt-in 3-D cuFFT path
and captures available NVIDIA libraries plus CUDA-aware MPI hints:

```
cd tests/profile/si64
NSTEPS=1 ./run_gpu_exploration.sh
```

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
  `WAVES_DEDPRO` MATMUL calls; set to `0` for the previous per-call copy path.
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
  `WAVES_ADDOPSI` on eligible non-stress paths. Non-superwave paths can keep
  OPSI resident from build and mass scaling; superwave paths currently build and
  mass-scale OPSI on the host before entering the resident region. The
  compatibility alias is `CPPAW_CUBLAS_ACC_OPSI_RESIDENCY`.
- `CPPAW_GPU_PSIM_PROPAGATE`: disabled by default. Set to `1` to run
  `WAVES$PROPAGATE` on the GPU for non-stress steps and copy the updated `PSIM`
  back before orthogonalization. The compatibility aliases are
  `CPPAW_GPU_PSIM_RESIDENCY` and `CPPAW_CUBLAS_ACC_PSIM_RESIDENCY`; they do not
  imply true cross-orthogonalization residency in the current implementation.
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
`cusolver_off`, one-rank CPU/NVHPC and eight-rank CPU/NVHPC references.

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
benchmark table is `combined_benchmark.tsv`.

Set `RUN_NVLAMATH=yes` to add a short NVLAMATH comparison suite,
`RUN_CUFFTW=yes` to add a short cuFFTW comparison suite, `RUN_CUFFT=yes` to add
a short native cuFFT comparison, `RUN_GPU_ACC=yes` to compare one-rank GPU
against eight-rank CPU/OpenBLAS and CPU/NVHPC references, or
`RUN_CUSOLVER=yes` to add the same one-rank cuSOLVER versus eight-rank CPU/NVHPC
resource comparison.

The overnight defaults use the recommended production-style cases:
`MAIN_CASES="nvhpc_cpu cublas cublas_off"`,
`SCALING_CASES="nvhpc_cpu cublas"`,
`GPU_ACC_CASES="cpu nvhpc_cpu gpu_resident gpu_resident_nosync gpu gpu_off"`
and `THRESHOLDS="1e7"`. Set
`RUN_GPU_DIAGNOSTICS=yes` to add `gpu_nosync`, `gpu_force_all` and the
`gpu_no_*` ablation cases, set `RUN_BAND_BENCHMARK=yes` to add the
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
