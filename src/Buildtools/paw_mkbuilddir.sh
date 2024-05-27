#!/bin/bash
################################################################################
##  paw_mkbuilddir fills the Build directory with symbolic links to
##  the source files. The modification times of the links is set identical 
##  to that of the source files.
################################################################################
#-------------------------------------------------------------------------------
#  help message
#-------------------------------------------------------------------------------
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \t paw_mkbbuilddir.sh options \n"
USAGE="$USAGE installation script for the cppaw package\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -b cppaw base directory\n"
USAGE="$USAGE \t -o build directory\n"
USAGE="$USAGE \t -i list of include files in a single string\n"
USAGE="$USAGE \t -v verbose (false)\n"
USAGE="$USAGE \t -h prints this help message\n"
USAGE="$USAGE \n"

#-------------------------------------------------------------------------------
#  Resolve arguments
#-------------------------------------------------------------------------------
export VERBOSE=false
export BUILDDIR=
export BASEDIR=
export INCLUDES=
while getopts :b:i:o:vh OPT ; do
  case $OPT in
    b) BASEDIR=$OPTARG ;;
    i) INCLUDES=$OPTARG ;;
    o) BUILDDIR=$OPTARG ;;
    v) VERBOSE=true ;;
    h) echo -e $USAGE ; exit 0  ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0" >&2
      echo "invalid option -$OPTARG" >&2
      echo "retrieve argument list with:" >&2
      echo "$0 -h" >&2
      exit 1
      ;;
    :)    # no argument passed
      ;;
  esac
done

if [[ -z ${BASEDIR} ]] ; then
  echo "error in $0: BASEDIR not specified (option -i missing)"
  exit 1
else
  if [[ ! -d ${BASEDIR}/src ]] ; then
    echo "error in $0: BASEDIR specified with option -i is invalid"
    echo "no subdirectory src"
    echo "specified BASEDIR=${BASEDIR}"
    exit 1
  fi 
fi
if [[ -z ${BUILDDIR} ]] ; then
  echo "error in $0: BUILDDIR not specified (option -o missing)"
  exit 1
fi

#-------------------------------------------------------------------------------
#-- cread Build directory
#-------------------------------------------------------------------------------
if [[ ! -d ${BUILDDIR} ]]     ; then mkdir ${BUILDDIR} ; fi
if [[ ! -d ${BUILDDIR}/etc ]] ; then mkdir ${BUILDDIR}/etc ; fi
if [[ ! -d ${BUILDDIR}/doc ]] ; then mkdir ${BUILDDIR}/doc ; fi

#-------------------------------------------------------------------------------
#-- include files
#-------------------------------------------------------------------------------
for X in ${INCLUDES}; do
   SOURCE=${X}
   TARGET=${BUILDDIR}/${X##*/}
   ln -sf ${SOURCE} ${TARGET}
   touch -hr ${SOURCE} ${TARGET}
done

#-------------------------------------------------------------------------------
#-- administration codes
#-------------------------------------------------------------------------------
export ADMINLIST=""
ADMINLIST="${ADMINLIST} Buildtools/dollar_ok.sh"
ADMINLIST="${ADMINLIST} Buildtools/f90pp.in"
ADMINLIST="${ADMINLIST} Buildtools/f90pp.sed"
ADMINLIST="${ADMINLIST} Buildtools/f90pp_tmplts.f90"
ADMINLIST="${ADMINLIST} Buildtools/parmfilewriter.f90"
ADMINLIST="${ADMINLIST} Buildtools/paw_versioninfo.sh"

for X in ${ADMINLIST}; do
   SOURCE=${BASEDIR}/src/${X}
   TARGET=${BUILDDIR}/etc/${X##*/}
   ln -sf ${SOURCE} ${TARGET}
   touch -hr ${SOURCE} ${TARGET}
done

#--- cppaw_version.info for releases -------------------------------------------
X=${BASEDIR}/cppaw_version.info
if [[ -e ${X} ]] ; then
   SOURCE=${X}
   TARGET=${BUILDDIR}/etc/${X##*/}
   ln -sf ${SOURCE} ${TARGET}
   touch -hr ${SOURCE} ${TARGET}
fi

if [[ ${VERBOSE} = true ]] ; then
   echo ADMINLIST=${ADMINLIST}
   ls -Ghal ${BUILDDIR}/etc/
fi


#-------------------------------------------------------------------------------
#-- object files derived from FLIBLIST and  will become part                  --
#-- of the paw library libpaw.a.                                              --
#-------------------------------------------------------------------------------
export LIBLIST=" \
        paw_trace \
        paw_error \
        paw_filehandler \
        paw_clock \
        paw_lock \
        paw_timing \
        paw_linkedlist \
        paw_strings \
        paw_dft \
        paw_dftaddendum \
        paw_constants \
        paw_spherical \
        paw_generalpurpose \
        paw_report \
        paw_periodictable \
        paw_radial \
        paw_schroedinger \
        paw_atomlib \
        paw_specialfunctions \
        paw_usage \
        paw_selftest \
        paw_strcio \
        paw_cell \
        paw_pdos \
        paw_banddata \
        paw_library \
        paw_polynom \
        paw_dimer  \
        paw_lmtobasics \
        paw_debug \
        paw_brillouin \
        paw_gaussian \
        paw_mpelib \
        paw_version "

#-------------------------------------------------------------------------------
#--  objects specific for the simulation code                                 --
#-------------------------------------------------------------------------------
export PAWLIST="paw_driver \
               paw_thermostat \
               paw_isolate \
               paw_assist \
               paw_lists \
               paw_constraints \
        paw_lmto \
        paw_simplelmto \
        paw_dmft \
        paw_fft \
        paw_augmentation \
        paw_softcore \
        paw_classical \
        paw_forcefield \
        paw_atoms \
        paw \
        paw_efg \
        paw_ioroutines \
        paw_iotra \
        paw_ionew \
        paw_qmmm \
        paw_cosmo \
        paw_warmup \
        paw_optfric \
        paw_setups \
        paw_potential \
        paw_occupations \
        paw_pairpotential \
        paw_graphics \
        paw_waves1 \
        paw_waves2 \
        paw_mixer \
        paw_cg \
        paw_kpoints \
        paw_vext \
        paw_vdw \
        paw_ci \
        paw_opteels"


# DIRTOOLS contains the names of the source files for the tools including 
# the path relative to the tool directory $(BASEDIR)/src/
DIRTOOLS="Tools/FromPOSCAR/paw_fromposcar \
         Tools/Grab/paw_grab \
         Tools/Murnaghan/paw_murnaghan \
         Tools/PDoS/paw_dos \
         Tools/PDoS/paw_dosplot \
         Tools/Polyhedra/paw_polyhedra \
         Tools/Preopt/paw_preopt \
         Tools/Stpa/paw_stpa \
         Tools/Stpa/paw_stpreport \
         Tools/Strc/paw_strc \
         Tools/Strc/paw_tostrc \
         Tools/Strc/paw_toxyz \
         Tools/Tra/paw_cleantra \
         Tools/Tra/paw_converttra \
         Tools/Tra/paw_tra \
         Tools/Wave/paw_1davpot \
         Tools/Wave/paw_cmcwave \
         Tools/Wave/paw_wave \
         Tools/Bands/paw_bands"

#-------------------------------------------------------------------------------
#-- fortran codes
#-------------------------------------------------------------------------------
export LIST=
LIST="${LIST} ${LIBLIST}"
LIST="${LIST} ${PAWLIST}"
LIST="${LIST} ${DIRTOOLS}"

for X in $LIST; do
   SOURCE=${BASEDIR}/src/${X}.f90
   TARGET=${BUILDDIR}/${X##*/}.f90pp
   ln -sf ${SOURCE} ${TARGET}
   touch -hr ${SOURCE} ${TARGET}
done

X=slatec
SOURCE=${BASEDIR}/src/${X}.f
TARGET=${BUILDDIR}/${X##*/}.f
ln -sf ${SOURCE} ${TARGET}
touch -hr ${SOURCE} ${TARGET}


if [[ ${VERBOSE} = true ]] ; then
   echo LIST=${LIST}
   ls -Ghal ${BUILDDIR}/*.f90pp
fi

#-------------------------------------------------------------------------------
#-- Scripts
#-------------------------------------------------------------------------------
export SCRIPTS=${BASEDIR}/src/Tools/Scripts/*.sh
for X in ${SCRIPTS}; do
   SOURCE=${X}
   TARGET=${BUILDDIR}/${X##*/}
   ln -sf ${SOURCE} ${TARGET}
   touch -hr ${SOURCE} ${TARGET}
done

if [[ ${VERBOSE} = true ]] ; then
   echo SCRIPTS=${SCRIPTS##*/}
  ls -Ghal ${BUILDDIR}/*.sh
fi
#-------------------------------------------------------------------------------
#-- Documentation
#-------------------------------------------------------------------------------
export DOCLIST=
DOCLIST="${DOCLIST} Docs/manual.tex"
DOCLIST="${DOCLIST} Docs/Figs"
DOCLIST="${DOCLIST} Docs/all.bib"
DOCLIST="${DOCLIST} Docs/doc.bib"
for X in ${DOCLIST}; do
  SOURCE=${BASEDIR}/src/${X}
  TARGET=${BUILDDIR}/doc/${X##*/}
  ln -sf ${SOURCE} ${TARGET}
  touch -hr ${SOURCE} ${TARGET}
done

if [[ ${VERBOSE} = true ]] ; then
  echo DOCLIST=${DOCLIST}
  ls -Ghal ${BUILDDIR}/doc
fi
