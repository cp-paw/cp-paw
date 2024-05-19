#!/bin/bash
export THISDIR=$(pwd)
#________1_________2_________3_________4_________5_________6_________7_________8
################################################################################
##     
##    paw_build.sh   
##     
##    build script for the cp-paw package 
##     
##     
##   Author P. Bloechl May, 2024 
################################################################################
#
#-------------------------------------------------------------------------------
#  help message
#-------------------------------------------------------------------------------
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \t paw_build.sh options \n"
USAGE="$USAGE installation script for the cppaw package\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -c choice\n"
USAGE="$USAGE \t -f parmfile \n"
USAGE="$USAGE \t -s suffix \n"
USAGE="$USAGE \t -v verbose (false)\n"
USAGE="$USAGE \t -h prints this help message\n"
USAGE="$USAGE \n"
USAGE="$USAGE -c is the suffix if no parmfile is present or -f ''\n"

#-------------------------------------------------------------------------------
#  Resolve arguments
#-------------------------------------------------------------------------------
export PARMFILE=
export SELECT=
export SUFFIX=
export VERBOSE=false
export PARALLEL=false
while getopts :c:f:s:vh OPT ; do
  case $OPT in
    c) SELECT=$OPTARG ;;
    f) PARMFILE=$OPTARG ;;
    s) SUFFIX=$OPTARG ;;
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

if [[ -z ${PARMFILE} ]] ; then
  PARMFILE=${THISDIR}/src/Buildtools/defaultparmfile
  echo "Warning from $0: no parmfile specified, setting default"
  echo "default parmfile= ${PARMFILE}"
fi

if [[ ! -f ${PARMFILE} ]] ; then
  echo "error in $0: parmfile does not exist"
  echo "             PARMFILE=${PARMFILE}"
  exit 1
fi

################################################################################
##
##  include $PARMFILE 
##
##  define: SUFFIX, PARALLEL
##          MAKE,AR,CPP,FC.LD
##          CPPFLAGS,FFLAGS,LDFLAGS
##          INCLUDES,LIBS,
##          BASEDIR,BUILDDIR,BINDIR
##
################################################################################
PROTECT_SUFFIX=$SUFFIX

source ${PARMFILE}

if [[ -n $PROTECT_SUFFIX ]] ; then
  echo "Warning: option -s to paw_build.sh overwrites value from parmfile"
  SUFFIX=$PROTECT_SUFFIX
fi

################################################################################
##  strip extra spaces
################################################################################
INCLUDES="$(echo ${INCLUDES} | tr -s '[:blank:]')"  # strip extra spaces
LIBS="$(echo ${LIBS} | tr -s '[:blank:]')"          # strip extra spaces
CPPFLAGS="$(echo ${CPPFLAGS} | tr -s '[:blank:]')"  # strip extra spaces
FFLAGS="$(echo ${FFLAGS} | tr -s '[:blank:]')"      # strip extra spaces
LDFLAGS="$(echo ${LDFLAGS} | tr -s '[:blank:]')"    # strip extra spaces

################################################################################
##  report settings
################################################################################
if [[ $VERBOSE = "true" ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report variables after parmfile---------------------------"
   echo "----------------------------------------------------------------------"
  LIST="SUFFIX PARALLEL MAKE AR CPP FC LD \
        CPPFLAGS FFLAGS LDFLAGS \
        INCLUDES LIBS \
        BASEDIR BUILDDIR BINDIR"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y="$(echo $Y | tr -s '[:blank:]')"  # strip extra spaces
      X="${X}.............................."
      X="${X:0:10}"
      echo -e "${X}: ${Y}" 
   done
   echo "-------------------------------------------------------done-----------"
fi

################################################################################
##     analyze system: MAKE, CPP, AR, FC, LD 
################################################################################
#-------------------------------------------------------------------------------
# operating system
#-------------------------------------------------------------------------------
# operating system may be Darwin (=macOS), Linux
OS="$(uname -s)"

#-------------------------------------------------------------------------------
# latexmk is required to produce the documentation
#-------------------------------------------------------------------------------
if [[ -z $(which latexmk) ]] ; then
  echo "error in $0: latexmk not installed or accessible"
  echo "see https://ctan.org/pkg/latexmk"
  echo "usually contained in the tex installation"
  echo "such as TeX Live https://tug.org/texlive/"
  exit 1
fi

#-------------------------------------------------------------------------------
# Make: must be gnu make, version 4.3 or later
# "grouped targets" were introduced only in version 4.3
#-------------------------------------------------------------------------------
#_____check if MAKE is defined__________________________________________________
if [[ -z ${MAKE} ]] ; then
  echo "error in $0: parameter MAKE not specified"
  exit 1
fi
#_____check if $MAKE exists and is executable___________________________________
if [[ ! -x ${MAKE} ]] ; then
  echo "error in $0: MAKE is not an executable file"
  exit 1
fi
#__ check if it is gnu make version 4.3 or later________________________________
VERSION=$(${MAKE} -v)
if [[ -z $(echo $VERSION | grep -Eo 'GNU Make') ]] ; then 
  echo "error in $0: MAKE is not GNU Make"
  exit 1
fi
# the version string contains two version numbers.
# outer parenthesis turns the result into a string array X[0] X[1] ...
X=($(echo ${VERSION} |  grep -Eo '[0-9]+\.[0-9]+'))
ID=${X[0]} 
MAJOR=${ID%.*}
MINOR=${ID#*.}
if [[ $MAJOR -lt 4 ]] ; then 
  echo "error in $0: MAKE is not version 4.3 or later"
  exit 1
fi
if [[ $MAJOR -eq 4 && $MINOR -lt 3 ]] ; then 
  echo "error in $0: MAKE is not version 4.3 or later"
  exit 1
fi

#-------------------------------------------------------------------------------
#   C - preprocessor
#-------------------------------------------------------------------------------
if [[ ! -x ${CPP} ]]; then
  echo "error in $0: no C-preprocessor found"
  echo 'specify variable CPP via parmfile'
  exit 1
fi

#-------------------------------------------------------------------------------
#   Archiver  (ar may become obsolete. switch to  tar instead?
#-------------------------------------------------------------------------------
if [[ ! -x ${AR} ]]; then
  echo "error in $0: no Archiver (ar) found"
  echo "specify variable AR via parmfile"
  exit 1
fi

#-------------------------------------------------------------------------------
# Fortran compiler
#-------------------------------------------------------------------------------
if [[ -z ${FC} ]] ; then
  echo "error in $0: parameter FC not specified"
  exit 1
fi
#_____check if $FC exists and is executable___________________________________
if [[ ! -x ${FC} ]] ; then
  echo "error in $0: FC is not an executable file"
  exit 1
fi

#-------------------------------------------------------------------------------
# Linker
#-------------------------------------------------------------------------------
if [[ -z ${LD} ]] ; then
  echo "error in $0: parameter LD not specified"
  exit 1
fi
#_____check if $FC exists and is executable___________________________________
if [[ ! -x ${LD} ]] ; then
  echo "error in $0: LD is not an executable file"
  exit 1
fi

#-------------------------------------------------------------------------------
# PARALLEL
#-------------------------------------------------------------------------------
if [[ ! ( $PARALLEL = true || $PARALLEL = false ) ]] ; then
  echo "error in $0: PARALLEL neither true nor false"
  exit 1
fi 

#-------------------------------------------------------------------------------
# CPPFLAGS
#-------------------------------------------------------------------------------
if [[ $PARALLEL = true && \
      -z "$(echo ${CPPFLAGS} | grep -Eo "CPPVARIABLE_PARALLEL")" ]] ; then
  CPPFLAGS="${CPPFLAGS} -DCPPVARIABLE_PARALLEL"
  echo "Warning: cpp flag "-DCPPVARIABLE_PARALLEL" required in parallel mode."
  echo "         Flag has been added to CPPFLAGS"
fi

if [[ -z $(echo ${CPPFLAGS} | grep -Eo "CPPVAR_FFT_FFTW3") ]] ; then
  CPPFLAGS="${CPPFLAGS} -DCPPVAR_FFT_FFTW3"
  echo "CPPFLAGS= $CPPFLAGS"
  echo "Warning: cpp flag -DCPPVAR_FFT_FFTW3 is mandatory"
  echo "         Flag has been added to CPPFLAGS"
fi

#CPPVAR_FEAST
#CPPVAR_SLEPC requires <FINCLUDE.SLEPCVEPSDEF.H> and slepceps.mod 
#                                                and  ISO_C_BINDING.mod
#CPPVAR_JADAMILU
#PETSC_USE_FORTRAN_DATATYPES

#-------------------------------------------------------------------------------
# FCFLAGS
#-------------------------------------------------------------------------------
# no test. can be anything, even empty
#-------------------------------------------------------------------------------
# LDFLAGS
#-------------------------------------------------------------------------------
# no test. can be anything, even empty

#-------------------------------------------------------------------------------
# LIBS
#-------------------------------------------------------------------------------
if [[ -z $(echo ${LIBS} | grep -Eo "fftw3") ]] ; then
  echo "error in $0: no fftw library specified on LIBS"
  exit
fi
if [[ -z $(echo ${LIBS} | grep -Eo "mpi") ]] ; then
  if [[ PARALLEL = true ]] ; then
    echo "error in $0: no mpi library specified in parallel mode"
    exit
  fi
fi
if [[ -z $(echo ${LIBS} | grep -Eo "xcf03") ]] ; then
  if [[ -z $(echo ${CPPFLAGS} | grep "NOLIBXC") ]] ; then
    echo "error in $0: no LIBXC library specified, while not switched off"
    exit
  fi
fi

for X in ${LIBS} ; do
  if [[ ${X:0:2} = -L ]] ; then
    DIR="${X#-L}"
  elif [[ ${X:0:2} = -l ]] ; then
    NAME="${X#-l}"
    LIB=${DIR}/lib${NAME}*
    if [[ -z  $LIB ]] ; then
      echo "error in $0: Libary $LIB doesnot exist"
      exit
    fi
  else
    echo "error in $0: Synatx error on LIBS"
    echo "all entries must have prefix -L or -l not separated by a blank"
  fi
done

#-------------------------------------------------------------------------------
# INCLUDES
#-------------------------------------------------------------------------------
if [[ -z $(echo ${INCLUDES} | grep -Eo "fftw3.f03") ]] ; then
  echo "error in $0: no fftw.f03 include file on INCLUDES"
  exit 1
fi

if [[ -z $(echo ${INCLUDES} | grep -Eo "xc_f03_lib_m.mod") ]] ; then
  echo "error in $0: no xc_f03_lib_m.mod module file on INCLUDES"
  exit 1
fi

if [[ -z $(echo ${INCLUDES} | grep -Eo "mpi_f08.mod") ]] ; then
  if [[ $PARALLEL = true ]] ; then
    echo "error in $0: no mpi_f08.mod module file on INCLUDES"
    exit 1
  fi
fi
for X in $INCLUDES ; do
  if [[ ! -f $X ]] ; then
    echo "error in $0: file $X in INCLUDES does not exist"
    echo "INCLUDES=$INCLUDES"
    exit 1
  fi 
done

################################################################################
##     Naming and locations                                                   ##
################################################################################
if [[ -z ${BASEDIR} ]] ; then
    echo "error in $0: BASEDIR not defined"
    exit
fi

if [[ ! -d ${BASEDIR} ]] ; then
  echo "error in $0: BASEDIR  does not exist"
  exit
fi
for X in src src/Buildtools src/Tools src/Docs ; do
  if [[ ! -d ${BASEDIR}/$X ]] ; then
    echo "error in $0: BASEDIR/$X  does not exist"
    exit
  fi
done
#__________ BUILDDIR____________________________________________________________
if [[ -z ${BUILDDIR} ]] ; then
  echo "error in $0: BUILDIR not defined"
  exit
fi
if [[ ! -d ${BUILDDIR} ]] ; then
  mkdir -p ${BUILDDIR} ; RC=$?
  if [[ $RC -ne 0 ]] ; then
    echo "error in $0: could not create BUILDIR"
    echo "BULDDIR=$BULDDIR"
    exit 1
  fi
fi
for X in etc doc ; do
   if [[ ! -d ${BUILDDIR}/$X ]] ; then
     mkdir ${BUILDDIR}/$X ; RC=$?
     if [[ $RC -ne 0 ]] ; then
        echo "error in $0: could not create BUILDIR/$X"
        echo "BULDDIR=$BULDDIR"
        exit 1
     fi
   fi
done

#________BINDIR_________________________________________________________________
if [[ -z ${BINDIR} ]] ; then
  echo "error in $0: BUILDIR not defined"
  exit
fi
if [[ ! -d ${BINDIR} ]] ; then
  mkdir ${BINDIR} ; RC=$?
  if [[ $RC -ne 0 ]] ; then
    echo "error in $0: could not create BINDIR"
    echo "BINDIR=$BINDIR"
    exit 1
  fi
fi
for X in include ; do
  if [[ ! -d ${BINDIR}/$X ]] ; then
    mkdir ${BINDIR}/$X ; RC=$?
    if [[ $RC -ne 0 ]] ; then
      echo "error in $0: could not create BINDIR/$X"
      echo "BINDIR/lib=${BINDIR}"
      exit 1
    fi
  fi
done

#________DOCDIR____(installation directory for the manal.pdf____________________
if [[ -z ${DOCDIR} ]] ; then
  echo "error in $0: DOCDIR not defined"
  exit
fi
if [[ ! -d ${DOCDIR} ]] ; then
  mkdir ${DOCDIR} ; RC=$?
  if [[ $RC -ne 0 ]] ; then
    echo "error in $0: could not create DOCDIR"
    echo "DOCDIR=$DOCDIR"
    exit 1
  fi
fi

################################################################################
##  strip extra spaces
################################################################################
INCLUDES=$(echo ${INCLUDES} | tr -s '[:blank:]')  # strip extra spaces
LIBS=$(echo ${LIBS} | tr -s '[:blank:]')  # strip extra spaces
CPPFLAGS=$(echo ${CPPFLAGS} | tr -s '[:blank:]')  # strip extra spaces
FFLAGS=$(echo ${FFLAGS} | tr -s '[:blank:]')  # strip extra spaces
LDFLAGS=$(echo ${LDFLAGS} | tr -s '[:blank:]')  # strip extra spaces

################################################################################
##     write parms.in_use
################################################################################
export LIST="SUFFIX PARALLEL\
             MAKE AR CPP FC LD\
             CPPFLAGS FFLAGS LDFLAGS \
             LIBS INCLUDES\
             BASEDIR BUILDDIR BINDIR"
#___write parameters to a temporary file $TMP and copy only if it differs_______
#___from ${BUILDDIR}/etc/parms.in_use___________________________________________
TMP=$(mktemp)
for X in $LIST ; do
  eval "Y=\${$X}"   # Y=value of the variable with name $X
  echo "export $X=$Y" >> $TMP
done
#___copy parms.in_use only if parameters have changed___________________________
#___to avoid a Makefile cascade_________________________________________________
diff "$TMP" "${BUILDDIR}/etc/parms.in_use" 
RC=$?    # RC=0: files are identical, RC=1: files differ; RC>1: diff fails
         # RC differs from zero also, when a file is missing
if [[ $RC -ne 0 ]] ; then
  cp ${TMP} ${BUILDDIR}/etc/parms.in_use
fi
rm -f $TMP

################################################################################
##     construct documentation
################################################################################
export MOPTS="-j 10"
#export MOPTS=""

export PARMLIST="DOCDIR BASEDIR"
export SEDCOMMANDS=$(mktemp)
for X in ${PARMLIST}; do
  eval "Y=\${$X}"
  echo "s|@${X}@|${Y}|g" >> ${SEDCOMMANDS}
done
#_____construct Makefile in build directory_____________________________________
sed -f $SEDCOMMANDS ${BASEDIR}/src/Buildtools/makedocs.in > ${BUILDDIR}/doc/Makefile
rm -f ${SEDCOMMANDS}

echo "................................................................make docs"
(cd ${BUILDDIR}/doc &&  ${MAKE} ${MOPTS})
echo "................................................................made docs"

################################################################################
##     compile
################################################################################
#
#______________________write sed file___________________________________________
export PARMLIST="SUFFIX BASEDIR \
                 MAKE AR FC LD CPP \
                 CPPFLAGS FFLAGS LDFLAGS \
                 LIBS INCLUDES"
export SEDCOMMANDS=$(mktemp)
for X in ${PARMLIST}; do
  eval "Y=\${$X}"
  echo "s|@${X}@|${Y}|g" >> ${SEDCOMMANDS}
done
#_____construct Makefile in build directory_____________________________________
sed -f $SEDCOMMANDS ${BASEDIR}/src/Buildtools/Makefile.in > ${BUILDDIR}/Makefile
rm -f ${SEDCOMMANDS}

#emacs ${BUILDDIR}/Makefile

# two calls to make are required because the first call constructs an 
# intermediate make file that considers module files
echo ".............................................................make prepare"
(cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} prepare)
RC=$?
if [[ $RC -ne 0 ]] ; then 
  echo "error in $0: make prepare in BUILDDIR exited with RC=$RC"
  exit 1
fi
echo ".............................................................made prepare"

if [[ ${PARALLEL} = true ]] ; then
  echo "....................................................... make executable"
  (cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} executable)
  RC=$?
  if [[ $RC -ne 0 ]] ; then 
    echo "error in $0: make executable in BUILDDIR exited with RC=$RC"
    exit 1
  fi
  echo " .......................................................made executable"
elif [[ ${PARALLEL} = false ]] ; then
  echo "...............................................................make all"
  (cd ${BUILDDIR} &&  ${MAKE} ${MOPTS}  all)
  RC=$?
  if [[ $RC -ne 0 ]] ; then 
    echo "error in $0: make all in BUILDDIR exited with RC=$RC"
    exit 1
  fi
  echo " ..............................................................made all"
fi

################################################################################
##     install
################################################################################
export PARMLIST="BINDIR SUFFIX PARALLEL"
export SEDCOMMANDS=$(mktemp)
for X in ${PARMLIST}; do
  eval "Y=\${$X}"
  echo "s|@${X}@|${Y}|g" >> ${SEDCOMMANDS}
done
#_____construct Makefile in build directory_____________________________________
sed -f $SEDCOMMANDS ${BASEDIR}/src/Buildtools/makeinstall.in > ${BUILDDIR}/install.mk
rm -f ${SEDCOMMANDS}

echo "............................................................make install"
if [[ $PARALLEL = false ]] ; then
  (cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} -f install.mk all)
  RC=$?
  if [[ $RC -ne 0 ]] ; then 
    echo "error in $0: make all in BUILDDIR exited with RC=$RC"
    exit 1
  fi
elif [[ $PARALLEL = true ]] ; then
  (cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} -f install.mk executable)
  RC=$?
  if [[ $RC -ne 0 ]] ; then 
    echo "error in $0: make all in BUILDDIR exited with RC=$RC"
    exit 1
  fi
fi
echo "................................................installation finished"
echo "cppaw installed in ${BINDIR}:"
# echo "files in ${BINDIR}:"
# ls ${BINDIR}
# echo "files in ${BINDIR}/include:"
# ls ${BINDIR}/include

