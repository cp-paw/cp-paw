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

- fortran compiler (e.g. gfortran, ifort)
- pkg-config, latexmk, GNU Make (version 4.3 or later)
- bash, cpp, ar
- tex (latex) distribution (e.g. TeX Live) 
- LAPACK, BLAS, FFTW3, MPI (optional), LIBXC (optional)
- tools: xmgrace, gnuplot, avogadro1

## Installation

1. Download the cppaw distribution from [Github](https://github.com/cp-paw/cp-paw).
2. Unpack the distribution in a directory. I will refer to this directory as the base directory.
3. Add the following lines to your profile (e.g. ~/.zshrc, ~/.bashrc, ~/.profile). Replace the definition of `PAWDIR` by the name of the base directory.
   ```
   export PAWDIR="path to your cppaw distribution"
   export PATH=${PAWDIR}/bin/fast:${PAWDIR}/bin/fast_parallel:{PAWDIR}/bin/dbg:${PATH}
   ```
4. Execute the following commands in the base directory:
   ```
   ./paw_install
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

## Optional Skala functional

The `skala` branch adds Skala 1.1 native-grid exchange-correlation evaluation
with joint periodic PAW reconstruction. It is based on `main` and does not
require the NVIDIA extensions in `cp-paw-nvhpc`. Conventional `dbg`, `fast`
and `fast_parallel` builds do not link PyTorch, FTorch or CUDA. Skala is used
only when explicitly selected in the input.

Provide a CPU PyTorch or LibTorch installation, then prepare the optional
Fortran bridge and hash-pinned model:

```sh
FC=gfortran src/Buildtools/paw_skala_setup.sh --device cpu --download-model
src/Buildtools/paw_build.sh -c skala_cpu_fast -j16 -z
src/Buildtools/paw_build.sh -c skala_cpu_fast_parallel -j16 -z
```

Set `CPPAW_SKALA_FTORCH_ROOT` for a custom bridge prefix. The installer builds
the additional CPU Skala binaries with `CPPAW_INSTALL_SKALA_CPU=require`.
The default is `no`.

Model inference can also use CUDA through a CUDA-enabled FTorch bridge,
independently of CP-PAW's NVIDIA linear-algebra and FFT backends. Optional
OpenACC grid and source kernels are retained, but require an OpenACC-capable
build. NVIDIA build profiles and conventional wave-function residency belong
to the separate `cp-paw-nvhpc` branch.

For model-free reconstruction and derivative tests after a normal build:

```sh
bash tests/unittests/skala_reconstruction/run.sh bin/Build_fast
```

See [SkalaBridge](src/SkalaBridge/README.md) for input settings, CPU and CUDA
installation, quadrature and cache controls, and optional Si2 and crystal
verification. Finite-difference checks establish derivative consistency,
not agreement of the functional implementation with an independent reference.

## Optional native SCAN functional

The native Fortran SCAN kernel can use the same joint periodic PAW
reconstruction and adjoint as Skala. It requires neither a model file nor
PyTorch, FTorch, LibXC or a GPU at runtime. A conventional `fast` or
`fast_parallel` build includes this opt-in path.

```text
!DFT TYPE=10
 !SCAN RADIALPOINTS=200 LEBEDEVEXACTNESS=53 CHECK=T !END
!END
```

Here `TYPE=10` selects the PBE reference for atomic setup construction.
The `SCAN` block replaces the solid or molecule XC energy and its complete
operator, including the positive kinetic-energy-density contribution.
It does not add SCAN on top of PBE. Collinear spin-polarized calculations
are supported. Native SCAN atomic setup generation is not provided.
The native kernel was contributed by Konstantin Tamoev and implements
[SCAN exchange and correlation](https://doi.org/10.1103/PhysRevLett.115.036402).
The verification drivers use the setup definitions in the supplied structure,
including built-in IDs or complete `!SPECIES!AUGMENT` blocks. They do not attach
the historical `PARMS_STP` file, which did not override the built-in setups.

The physical spin densities, gradients and positive kinetic energy densities
are reconstructed before the local functional is evaluated. Quadrature
weights multiply the resulting energy and derivatives, not its inputs.
The isolated frozen-core SCAN XC energy is subtracted as a constant to
retain CP-PAW's valence-energy convention. The protocol reports both the
reconstructed energy and this reference. Skala's existing energy convention
is unchanged.

This path integrates the functional of the jointly reconstructed fields.
It does not evaluate separate smooth, all-electron one-center and pseudo
one-center XC energies. Those two PAW approximations need not coincide at
finite partial-wave completeness. The joint path includes nonlinear cross
terms through the reconstructed fields and differentiates the same discrete
energy through the smooth fields, on-site density matrices and quadrature
geometry. It uses only first functional derivatives, not the conventional
one-center XC Taylor expansion.

The complete smooth, tau and one-center operators are mandatory. Native
SCAN rejects model/device settings, disabled operator contributions and
`VDW=T`, which would otherwise select the auxiliary PBE dispersion parameters.
Conventional input files without a `SCAN` or `SKALA` block are unchanged.

Run the kernel, independent LibXC-reference and backend tests with

```sh
bash tests/unittests/scan/run.sh bin/Build_fast
python3 tests/unittests/scan/input_checks.py bin/fast/paw_fast.x
bash tests/unittests/skala_reconstruction/run.sh bin/Build_fast
```

The reconstruction suite includes independent pointwise kinetic-density
references for complex Bloch waves, packed real orbitals, nonuniform
occupations, analytic core orbitals and complex s/p/d PAW superpositions.
These checks test the fields themselves, not just energy derivatives.

`tests/fulltests/scan_si2/run.py` audits an existing restart, while
`tests/fulltests/scan_si2/relax.py` requires explicit orbital-residual and
occupation-commutator convergence. The shared orbital, force and stress
verification drivers accept `--functional scan` without a model file.
Their finite-difference tests establish derivative consistency, not
independent physical accuracy. Joint reconstructed-density SCAN remains
under validation, particularly for quadrature and basis convergence and
comparison with independently converged calculations. Reference
comparisons must also account for the PAW dataset, frozen-core reference
and relativistic treatment. The retained
numerical second- and third-derivative kernel entry points are diagnostic
and are not used by this PAW path.

Start density-grid convergence checks at `CDUAL=4`, which is also the
default in the native SCAN verification drivers. `CHECK=T` reports low-density,
negative-tau and spin-channel von Weizsaecker-bound violations in the raw
reconstructed fields, including the quadrature-weighted electron count in
the affected region. These diagnostics can expose underresolved fields.
Increasing the atom-grid size alone does not refine the underlying density grid.
Integrated Pauli and negative-tau deficits report the magnitude of the
violations in Hartree using absolute quadrature weights. They distinguish
rounding-level sign changes from finite deficits. Neither affected electron
counts nor integrated deficits are error bounds for the XC energy.
The native kernel retains the smooth SCAN continuation at negative iso-orbital
indicator alpha rather than imposing an extra alpha=0 clamp. Its boundary
references use LibXC with Fermi-hole-curvature (FHC) enforcement disabled.
LibXC's optional FHC enforcement changes gradient invariants and is a distinct
regularization. Match this convention in independent reference calculations.

## License

The CP-PAW code is distributed under the GNU Public License Version 3.
See the LICENSE file.


