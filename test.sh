#!/bin/bash
PAWDIR=/Users/ptpb/Tree/PAW/devel
BINFAST=${PAWDIR}/bin/osx_gfortran_mkl/Objects/fast

cd ${BINFAST}

${PAWDIR}/bin/osx_gfortran_mkl/f90pp \
         -DCPPVAR_COMPILER_GFORTRAN -DCPPVAR_FFT_FFTW \
    ${PAWDIR}/src/paw_version.f90 \
    >${BINFAST}/paw_version_d1.f90

cp ${BINFAST}/paw_version_d1.f90 ${BINFAST}/paw_version_d.f90

${PAWDIR}/src/Buildtools/Version/getversion.sh ${BINFAST}/paw_version_d.f90
