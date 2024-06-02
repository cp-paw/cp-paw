CP-PAW code package
====================

See https://cppaw.org for further information (Currently, the
description on https://cppaw.org refers to an older release and does
not apply to the present implementation.)


Configuration and Installation instructions
===========================================

1) Preconditions

- fortran compiler (e.g. gfortran, ifort)
- pkg-config, latexmk, GNU Make, version 4.3 or later
- bash, cpp, ar, 
- A tex (latex) distribution (e.g. TeX Live) 
- libraries: LAPACK, BLAS, FFTW3, MPI (optional), LIBXC (optional)
  
For the tools, we also use xmgrace, gnuplot, avogadro1

2) Installation

- Download the cppaw distribution from Github.  

- Unpack the distribution in a directory. I will refer to this
  directory as the BASE directory

- add the following lines to your profile (e.g. ~/.zshrc, ~/.bashrc,
  ~/.profile). Replace the definition of PAWdir by the name of the
  base direcory

  export PAWDIR="path to your cppaw distribution"
  export PATH=${PAWDIR}/bin/fast:${PAWDIR}/bin/fast_parallel:${PAWDIR}/bin/dbg:${PATH}

- execute

./paw_install

  in the base directory

- if this does not work, make a copy of the default parmfile

cp src/Buildtools/defaultparmfile ./parmfile

After editing the parmfile open paw\_install and add "-f parmfile" as
argumnent for paw_build.sh, such as

src/Buildtools/paw_build.sh -f parmfile -c fast 
src/Buildtools/paw_build.sh -f parmfile -c fast_parallel

- consult the manual doc/manual.pdf for further information.

Copyright
=========
The CP-PAW code is distributed under the GNU Public License Version 3.
See the LICENSE file.



