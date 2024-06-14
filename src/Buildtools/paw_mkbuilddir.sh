#!/bin/bash
################################################################################
##  paw_mkbuilddir fills the Build directory with symbolic links to
##  the source files. The modification times of the links is set identical 
##  to that of the source files.
################################################################################
# option -h of ln: If the target_file or target_dir is a symbolic
# link, do not follow it.  This is most useful with the -f option, to
# replace a symlink which may point to a directory.

# option -f of ln: If the target file already exists, then unlink it
# so that the link may occur.
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
  echo "error in $0: BASEDIR not specified (option -i missing)"              >&2
  exit 1
else
  if [[ ! -d ${BASEDIR}/src ]] ; then
    echo "error in $0: BASEDIR specified with option -i is invalid"          >&2
    echo "no subdirectory src"                                               >&2
    echo "specified BASEDIR=${BASEDIR}"                                      >&2
    exit 1
  fi 
fi
if [[ -z ${BUILDDIR} ]] ; then
  echo "error in $0: BUILDDIR not specified (option -o missing)"             >&2
  exit 1
fi

#-------------------------------------------------------------------------------
#-- create Build directory
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
export ADMINLIST=$(${BASEDIR}/src/Buildtools/paw_srclist.sh -a)
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
#-- fortran codes
#-------------------------------------------------------------------------------
export LIST="$(${BASEDIR}/src/Buildtools/paw_srclist.sh -l) \
             $(${BASEDIR}/src/Buildtools/paw_srclist.sh -p) \
             $(${BASEDIR}/src/Buildtools/paw_srclist.sh -t)"
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
  if [[ -d ${TARGET} ]] ; then rm -f ${TARGET} ; fi
  ln -sf ${SOURCE} ${TARGET}
  touch -hr ${SOURCE} ${TARGET}
done

if [[ ${VERBOSE} = true ]] ; then
  echo DOCLIST=${DOCLIST}
  ls -Ghal ${BUILDDIR}/doc
fi
