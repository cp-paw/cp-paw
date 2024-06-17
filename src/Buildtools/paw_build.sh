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
USAGE="$USAGE \t -j jobs (jobs=nr parallel jobs. default 10)  \n"
USAGE="$USAGE \t -z suppress documentation manual.pdf\n"
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
export JOBS=10
export NODOC=false
while getopts :c:f:s:j:zvh OPT ; do
  case $OPT in
    c) SELECT=$OPTARG ;;
    f) PARMFILE=$OPTARG ;;
    s) SUFFIX=$OPTARG ;;
    j) JOBS=$OPTARG 
       export MOPTS="-j ${JOBS}"
       #  MOPTS="${MOPTS} --silent"
       MOPTS="${MOPTS} --debug=b"
       ;;
    z) NODOC=true ;;
    v) VERBOSE=true ;;
    h) echo -e $USAGE ; exit 0  ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0"                                                     >&2
      echo "invalid option -$OPTARG"                                         >&2
      echo "retrieve argument list with:"                                    >&2
      echo "$0 -h"                                                           >&2
      exit 1
      ;;
    :)    # no argument passed
      ;;
  esac
done

#  specify default parmfile, when no parmfile is explicitly specified
if [[ -z ${PARMFILE} ]] ; then
  PARMFILE=${THISDIR}/src/Buildtools/defaultparmfile
  echo "Using default parmfile ${PARMFILE}"
fi

if [[ ! -f ${PARMFILE} ]] ; then
  echo "error in $0: parmfile does not exist"                                >&2
  echo "             PARMFILE=${PARMFILE}"                                   >&2
  exit 1
fi

################################################################################
##
##  include $PARMFILE 
##
##  define: SUFFIX, PARALLEL
##          MAKE,AR,CPP,FC.LD
##          CPPFLAGS,FCFLAGS,LDFLAGS
##          INCLUDES,LIBS,
##          BASEDIR,BUILDDIR,BINDIR
##
################################################################################
PROTECT_SUFFIX=$SUFFIX

source ${PARMFILE}
RC=$?
if [[ RC -ne 0 ]] ; then 
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: parmfile returned with error"                           >&2
  echo "parmfile=${PARMFILE}"                                                >&2
  exit 1
fi

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
FCFLAGS="$(echo ${FCFLAGS} | tr -s '[:blank:]')"    # strip extra spaces
LDFLAGS="$(echo ${LDFLAGS} | tr -s '[:blank:]')"    # strip extra spaces

################################################################################
##  report settings
################################################################################
if [[ $VERBOSE = "true" ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report parameters after parmfile--------------------------"
   echo "----------------------------------------------------------------------"
   LIST="SUFFIX PARALLEL MAKE AR CPP FC LD \
        CPPFLAGS FCFLAGS LDFLAGS \
        INCLUDES LIBS \
        BASEDIR BUILDDIR BINDIR"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y="$(echo $Y | tr -s '[:blank:]')"  # strip extra spaces
      X="${X}.............................."
      X="${X:0:10}"
      echo -e "${X}: ${Y}" 
   done
   echo "-------------------------------------------parameter report done------"
fi

################################################################################
##     analyze system: MAKE, CPP, AR, FC, LD 
################################################################################
#  collect error messages before terminating the installation process
ERROR=false

#-------------------------------------------------------------------------------
# operating system
#-------------------------------------------------------------------------------
# operating system may be Darwin (=macOS), Linux
OS="$(uname -s)"

#-------------------------------------------------------------------------------
# latexmk is required to produce the documentation
#-------------------------------------------------------------------------------
if [[ ${NODOC} = false ]] ; then
  if [[ -z $(which latexmk) ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: latexmk not installed or accessible"                  >&2
    echo "see https://ctan.org/pkg/latexmk"                                  >&2
    echo "usually contained in the tex installation"                         >&2
    echo "such as TeX Live https://tug.org/texlive/"                         >&2
    ERROR=true
  fi
fi

#-------------------------------------------------------------------------------
# Make: must be gnu make, version 4.3 or later
# "grouped targets" were introduced only in version 4.3
#-------------------------------------------------------------------------------
#_____check if MAKE is defined__________________________________________________
if [[ -z ${MAKE} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: parameter MAKE not specified"                           >&2
  ERROR=true
fi
#_____check if $MAKE exists and is executable___________________________________
if [[ ! -x ${MAKE} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: MAKE is not an executable file"                         >&2
  ERROR=true
fi
#__ check if it is gnu make version 4.3 or later________________________________
VERSION=$(${MAKE} -v)
if [[ -z $(echo $VERSION | grep -Eo 'GNU Make') ]] ; then 
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: MAKE is not GNU Make"                                   >&2
  ERROR=true
fi
# the version string contains two version numbers.
# outer parenthesis turns the result into a string array X[0] X[1] ...
X=($(echo ${VERSION} |  grep -Eo '[0-9]+\.[0-9]+'))
ID=${X[0]} 
MAJOR=${ID%.*}
MINOR=${ID#*.}
if [[ $MAJOR -lt 4 ]] ; then 
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: MAKE is not version 4.3 or later"                       >&2
  ERROR=true
fi
if [[ $MAJOR -eq 4 && $MINOR -lt 3 ]] ; then 
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: MAKE is not version 4.3 or later"                       >&2
  ERROR=true
fi

#-------------------------------------------------------------------------------
#   C - preprocessor
#-------------------------------------------------------------------------------
if [[ ! -x ${CPP} ]]; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: no C-preprocessor found"                                >&2
  echo 'specify variable CPP via parmfile'                                   >&2
  ERROR=true
fi

#-------------------------------------------------------------------------------
#   Archiver  (ar may become obsolete. switch to  tar instead?
#-------------------------------------------------------------------------------
if [[ ! -x ${AR} ]]; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: no Archiver (ar) found"                                 >&2
  echo "specify variable AR via parmfile"                                    >&2
  ERROR=true
fi

#-------------------------------------------------------------------------------
# Fortran compiler
#-------------------------------------------------------------------------------
if [[ -z ${FC} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: parameter FC not specified"                             >&2
  ERROR=true
fi
#_____check if $FC exists and is executable___________________________________
if [[ ! -x ${FC} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: FC is not an executable file"                           >&2
  echo "FC=${FC}"                                                            >&2
  ERROR=true
fi

#-------------------------------------------------------------------------------
# Linker
#-------------------------------------------------------------------------------
if [[ -z ${LD} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: parameter LD not specified"                             >&2
  ERROR=true
fi
#_____check if $FC exists and is executable___________________________________
if [[ ! -x ${LD} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: LD is not an executable file"                           >&2
  echo "LD=${LD}"                                                            >&2
  ERROR=true
fi

#-------------------------------------------------------------------------------
# PARALLEL
#-------------------------------------------------------------------------------
if [[ ! ( $PARALLEL = true || $PARALLEL = false ) ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: PARALLEL neither true nor false"                        >&2
  echo "PARALLEL=${PARALLEL}"                                                >&2
  ERROR=true
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
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: no fftw library specified on LIBS"                      >&2
  ERROR=true
fi
if [[ -z $(echo ${LIBS} | grep -Eo "mpi") ]] ; then
  if [[ PARALLEL = true ]] ; then
    if [[ -z $(echo ${FC} | grep -Eo "mpi") ]] ; then
      echo "----------------------------------------------------------------">&2
      echo "error in $0: no mpi library specified in parallel mode"          >&2
      echo "no mpi compiler wrapper used: FC=${FC}"                          >&2
      ERROR=true
    fi
  fi
fi
if [[ -z $(echo ${LIBS} | grep -Eo "xcf03") ]] ; then
  if [[ -z $(echo ${CPPFLAGS} | grep "NOLIBXC") ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: no LIBXC library specified, while not switched off"   >&2
    ERROR=true
  fi
fi

for X in ${LIBS} ; do
  if [[ ${X:0:2} = -L ]] ; then
    DIR="${X#-L}"
  elif [[ ${X:0:2} = -l ]] ; then
    NAME="${X#-l}"
    LIB=${DIR}/lib${NAME}*
    if [[ -z  $LIB ]] ; then
      echo "----------------------------------------------------------------">&2
      echo "error in $0: Libary $LIB does not exist"                         >&2
      ERROR=true
    fi
  else
    echo "------------------------------------------------------------------">&2
    echo "error in $0: Syntax error on LIBS"                                 >&2
    echo "all entries must have prefix -L or -l not separated by a blank"    >&2
    ERROR=true
  fi
done

#  check whether commandline tools are installed on MacOS
if [[ $OS = Darwin ]] ; then
  if [[ -n $(xcode-select -p | grep -oE error) ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: command-line-tools of Xcode not installed"            >&2
    echo "install command-line tools of Xcode"                               >&2
    ERROR=true
  fi
fi

#-------------------------------------------------------------------------------
# INCLUDES
#-------------------------------------------------------------------------------
for X in $INCLUDES ; do
  if [[ ! -f $X ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: file $X in INCLUDES does not exist"                   >&2
    echo "INCLUDES=$INCLUDES"                                                >&2
    ERROR=true
  fi 
done

################################################################################
##     Naming and locations                                                   ##
################################################################################
if [[ -z ${BASEDIR} ]] ; then
    echo "error in $0: BASEDIR not defined"
    ERROR=true
fi

if [[ ! -d ${BASEDIR} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: BASEDIR  does not exist"                                >&2
  ERROR=true
fi
for X in src src/Buildtools src/Tools src/Docs ; do
  if [[ ! -d ${BASEDIR}/$X ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: BASEDIR/$X  does not exist"                           >&2
    ERROR=true 
  fi
done
#__________ BUILDDIR____________________________________________________________
if [[ -z ${BUILDDIR} ]] ; then
  echo "-----------------------------------------------------------------------"
  echo "error in $0: BUILDIR not defined"
  ERROR=true
fi
if [[ ! -d ${BUILDDIR} ]] ; then
  mkdir -p ${BUILDDIR} ; RC=$?
  if [[ $RC -ne 0 ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: could not create BUILDIR"                             >&2
    echo "BULDDIR=$BULDDIR"                                                  >&2
    ERROR=true
  fi
fi
for X in etc doc ; do
   if [[ ! -d ${BUILDDIR}/$X ]] ; then
     mkdir ${BUILDDIR}/$X ; RC=$?
     if [[ $RC -ne 0 ]] ; then
        echo "--------------------------------------------------------------">&2
        echo "error in $0: could not create BUILDIR/$X"                      >&2
        echo "BULDDIR=$BULDDIR"                                              >&2
        ERROR=true
     fi
   fi
done

#________BINDIR_________________________________________________________________
if [[ -z ${BINDIR} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: BUILDIR not defined"                                    >&2
  ERROR=true
fi
if [[ ! -d ${BINDIR} ]] ; then
  mkdir ${BINDIR} ; RC=$?
  if [[ $RC -ne 0 ]] ; then
    echo "------------------------------------------------------------------">&2
    echo "error in $0: could not create BINDIR"                              >&2
    echo "BINDIR=$BINDIR"                                                    >&2
    ERROR=true
  fi
fi
for X in include ; do
  if [[ ! -d ${BINDIR}/$X ]] ; then
    mkdir ${BINDIR}/$X ; RC=$?
    if [[ $RC -ne 0 ]] ; then
      echo "----------------------------------------------------------------">&2
      echo "error in $0: could not create BINDIR/$X"                         >&2
      echo "BINDIR/lib=${BINDIR}"                                            >&2
      ERROR=true
    fi
  fi
done

#________DOCDIR____(installation directory for the manal.pdf____________________
if [[ -z ${DOCDIR} ]] ; then
  echo "--------------------------------------------------------------------">&2
  echo "error in $0: DOCDIR not defined"                                     >&2
  ERROR=true
fi

################################################################################
##  strip extra spaces
################################################################################
INCLUDES=$(echo ${INCLUDES} | tr -s '[:blank:]')  # strip extra spaces
LIBS=$(echo ${LIBS} | tr -s '[:blank:]')  # strip extra spaces
CPPFLAGS=$(echo ${CPPFLAGS} | tr -s '[:blank:]')  # strip extra spaces
FCFLAGS=$(echo ${FCFLAGS} | tr -s '[:blank:]')  # strip extra spaces
LDFLAGS=$(echo ${LDFLAGS} | tr -s '[:blank:]')  # strip extra spaces

export ERRMSG=""
ERRMSG="${ERRMSG}\n SUFFIX  \t=${SUFFIX}"
ERRMSG="${ERRMSG}\n PARALLEL\t=${PARALLEL}"
ERRMSG="${ERRMSG}\n MAKE   \t\t=${MAKE}"
ERRMSG="${ERRMSG}\n AR     \t\t=${AR}"
ERRMSG="${ERRMSG}\n CPP    \t\t=${CPP}"
ERRMSG="${ERRMSG}\n FC     \t\t=${FC}"
ERRMSG="${ERRMSG}\n LD     \t\t=${LD}"         
ERRMSG="${ERRMSG}\n CPPFLAGS\t=${CPPFLAGS}"  
ERRMSG="${ERRMSG}\n FCFLAGS \t=${FCFLAGS}"   
ERRMSG="${ERRMSG}\n LDFLAGS \t=${LDFLAGS}"   
ERRMSG="${ERRMSG}\n LIBS   \t\t=${LIBS}"      
ERRMSG="${ERRMSG}\n INCLUDES\t=$INCLUDES"    
ERRMSG="${ERRMSG}\n DOCDIR  \t=${DOCDIR}"    
ERRMSG="${ERRMSG}\n BINDIR  \t=${BINDIR}"    
ERRMSG="${ERRMSG}\n BUILDDIR\t=${BUILDDIR}"  
ERRMSG="${ERRMSG}\n BASEDIR \t=${BASEDIR}"   
if [[ ${ERROR} != false ]] ; then
  echo "---------------report before error exit-----------------------------">&2
  echo -e ${ERRMSG}                                                          >&2
  echo "--------------------------------------------------------------------">&2
  echo "error exit from $0"                                                  >&2
  exit 1
fi

################################################################################
##     fill build directory
################################################################################
${BASEDIR}/src/Buildtools/paw_mkbuilddir.sh -b ${BASEDIR} -o${BUILDDIR} -i "$INCLUDES"
RC=$?
if [[ ${RC} -ne 0 ]] ; then 
  echo "error in $0: paw_mkbuilddir.sh returned with error code ${RC}"       >&2
  exit 1
fi

################################################################################
##     write parms.in_use to ${BUILDDIR}/etc/
################################################################################
export LIST="SUFFIX PARALLEL\
             MAKE AR CPP FC LD\
             CPPFLAGS FCFLAGS LDFLAGS \
             LIBS INCLUDES"
#___write parameters to a temporary file $TMP and copy only if it differs_______
#___from ${BUILDDIR}/etc/parms.in_use___________________________________________
TMP=$(mktemp)
echo "#!/bin/bash" >> $TMP
for X in $LIST ; do
  eval "Y=\${$X}"   # Y=value of the variable with name $X
  echo "export $X=\"$Y\"" >> $TMP
done
echo "export BASEDIR=./" >> $TMP
echo "export BINDIR=./bin/rebuild" >> $TMP
echo "export BUILDDIR=./bin/Build_rebuild" >> $TMP
echo "export DOCDIR=./doc" >> $TMP

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
if [[ ${DOCDIR}/manual.pdf -ot ${BASEDIR}/Docs/manual.tex ]] ; then
  NODOC=true
fi

if [[ ${NODOC} = false ]] ; then

  #__create doc directory_______________________________________________________
  if [[ ! -d ${DOCDIR} ]] ; then
    mkdir ${DOCDIR} ; RC=$?
    if [[ $RC -ne 0 ]] ; then
      echo "error in $0: could not create DOCDIR"                            >&2
      echo "DOCDIR=$DOCDIR"                                                  >&2
      echo "Potential racing condition of several paw_build.sh."             >&2
      echo "Consider using option -z of paw_build.sh \
            on the competing calls except one."                              >&2
      exit=1
    fi
  fi

  #__collect sed commands_______________________________________________________
  export PARMLIST="DOCDIR BASEDIR"
  export SEDCOMMANDS=$(mktemp)
  for X in ${PARMLIST}; do
    eval "Y=\${$X}"
    echo "s|@${X}@|${Y}|g" >> ${SEDCOMMANDS}
  done

  #_____construct Makefile in build directory___________________________________
  sed -f $SEDCOMMANDS ${BASEDIR}/src/Buildtools/makedocs.in \
                      > ${BUILDDIR}/doc/Makefile
  rm -f ${SEDCOMMANDS}

  echo "..............................................................make docs"
  (cd ${BUILDDIR}/doc &&  ${MAKE} ${MOPTS} )
  export RC=$?
  if [[ ${RC} -ne 0 ]] ; then
    echo "----------------------------------------------------------------" >&2
    echo "error in $0: making documentation failed"                         >&2
    echo "$(sed -n '/^[!l]/p' ${BUILDDIR}/doc/manual.log)"                  >&2
    echo "consult ${BUILDDIR}/doc/manual.log for details"                   >&2
    echo "shutting down $0"                                                 >&2
    echo "----------------------------------------------------------------" >&2
    exit 1
  fi
  echo "..............................................................made docs"
fi

################################################################################
##     compile
################################################################################

#  collect lists of source files (w/o extension, relative to src)
export LIBLIST=$(${BASEDIR}/src/Buildtools/paw_srclist.sh -l)
export ADMINLIST=$(${BASEDIR}/src/Buildtools/paw_srclist.sh -a)
export PAWLIST=$(${BASEDIR}/src/Buildtools/paw_srclist.sh -p)
export TOOLLIST=$(${BASEDIR}/src/Buildtools/paw_srclist.sh -t)
#
#______________________write sed file___________________________________________
export PARMLIST="SUFFIX BASEDIR \
                 MAKE AR FC LD CPP \
                 CPPFLAGS FCFLAGS LDFLAGS \
                 LIBS INCLUDES \
                 ADMINLIST LIBLIST PAWLIST TOOLLIST"
export SEDCOMMANDS=$(mktemp)
for X in ${PARMLIST}; do
  eval "Y=\${$X}"
  echo "s|@${X}@|${Y}|g" >> ${SEDCOMMANDS}
done
#_____construct Makefile in build directory_____________________________________
sed -f $SEDCOMMANDS ${BASEDIR}/src/Buildtools/Makefile.in > ${BUILDDIR}/Makefile
rm -f ${SEDCOMMANDS}

# two calls to make are required because the first call constructs an 
# intermediate make file that considers module files
echo ".............................................................make prepare"
(cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} prepare)
RC=$?
if [[ $RC -ne 0 ]] ; then 
  echo "error in $0: make prepare in BUILDDIR exited with RC=$RC"            >&2
  exit 1
fi
echo ".............................................................made prepare"

echo "..............................................................make big.mk"
(cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} -f big.mk)
export RC=$?
if [[ $RC -ne 0 ]] ; then 
  echo "error in $0: make prepare in BUILDDIR exited with RC=$RC"            >&2
  exit 1
fi
echo "..............................................................made big.mk"


echo ".................................................................make all"
(cd ${BUILDDIR} &&  ${MAKE} ${MOPTS}  all)
export RC=$?
if [[ $RC -ne 0 ]] ; then 
  echo "error in $0: make all in BUILDDIR exited with RC=$RC"                >&2
  exit 1
fi
echo " ................................................................made all"

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
(cd ${BUILDDIR} &&  ${MAKE} ${MOPTS} -f install.mk all)
RC=$?
if [[ $RC -ne 0 ]] ; then 
  echo "error in $0: make all in BUILDDIR exited with RC=$RC"                >&2
  exit 1
fi
echo ".............................................................made install"


if [[ ${PARALLEL} = true && $(uname -s) = Darwin ]] ; then
  echo "...........................................................code signing"
  ${BASEDIR}/src/Buildtools/Codesign/paw_codesign.sh ${BINDIR}/ppaw_${SUFFIX}.x
  echo "...........................................................code signed"
fi
echo "-------------------------------------------------------------------------"
echo "---------- cppaw installed in ${BINDIR} -----"
echo "-------------------------------------------------------------------------"

# echo "files in ${BINDIR}:"
# ls ${BINDIR}
# echo "files in ${BINDIR}/include:"
# ls ${BINDIR}/include

