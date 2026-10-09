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

## License

The CP-PAW code is distributed under the GNU Public License Version 3.
See the LICENSE file.


