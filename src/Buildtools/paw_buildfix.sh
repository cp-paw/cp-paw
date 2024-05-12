#!/bin/bash
export THISDIR=$(pwd)
#________1_________2_________3_________4_________5_________6_________7_________8
################################################################################
##     
##    paw_buildfix.sh   
##     
##    build script for the cp-paw package 
##     
##     
##   Author P. Bloechl May, 2024 
################################################################################

################################################################################
##  report settings before fix
################################################################################
if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report variables after parmfile---------------------------"
   echo "----------------------------------------------------------------------"
  LIST="SUFFIX PARALLEL MAKE AR CPP FC LD \
        CPPFLAGS FFLAGS LDFLAGS \
        INCLUDES LIBS \
        BASEDIR BUILDDIR BINDIR"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:20}
      echo -e "$X: $Y" 
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
OS=$(uname -s)
echo OS=$OS


#-------------------------------------------------------------------------------
# Make: must be gnu make, version 4.3 or later
# "grouped targets" were introduced only in version 4.3
#-------------------------------------------------------------------------------
for X in make gmake ; do
  # exit loop if MAKE has already been found
  if [[ -n ${MAKE} ]] ; then break ; fi
  MAKE=$(which $X) 
  if [[ -z $MAKE ]] ; then continue ; fi
  VERSION=$(${MAKE} -v)
  if [[ -z $(echo $VERSION | grep -Eo 'GNU Make') ]] ; then 
    MAKE=""; continue 
  fi
  # the version string contains two version numbers.
  # outer parenthesis turns the result into a string array X[0] X[1] ...
  X=($(echo $VERSION |  grep -Eo '[0-9]+\.[0-9]+'))
  ID=${X[0]} 
  MAJOR=${ID%.*}
  MINOR=${ID#*.}
  if [[ $MAJOR -lt 4 ]] ; then MAKE="" ; continue ;fi
  if [[ $MAJOR -eq 4 && $MINOR -lt 3 ]] ; then MAKE=""; continue;  fi
done
if [[ ! -x $MAKE ]]; then
  echo "error in $0: no GNU Make version >= 4.3 found"
  echo "MAKE=$MAKE"
  echo 'specify variable MAKE via parmfile'
  exit 1
fi

#-------------------------------------------------------------------------------
#   C - preprocessor
#-------------------------------------------------------------------------------
if [[ ! -x $CPP ]]; then
  CPP=$(which cpp) 
fi
if [[ ! -x $CPP ]]; then
  echo "error in $0: no C-preprocessor found"
  echo 'specify variable CPP via parmfile'
  exit 1
fi

#-------------------------------------------------------------------------------
#   Archiver  (ar may become obsolete. switch to  tar instead?
#-------------------------------------------------------------------------------
if [[ ! -x $AR ]] ; then
  AR=$(which ar)
fi
if [[ ! -x $AR ]]; then
  echo "error in $0: no Archiver (ar) found"
  echo "specify variable AR via parmfile"
  exit 1
fi

#-------------------------------------------------------------------------------
# Fortran compiler
#-------------------------------------------------------------------------------
for X in ifx xlf flang ifort gfortran ; do
  # exit loop if FC has already been found
  if [[ -x ${FC} ]] ; then break ; fi
  FC=$(which $X) 
  if [[ -z $FC ]] ; then continue ; fi
done
if [[ ! -x $FC ]] ; then 
  echo "error in $0: no fortran compiler found"
  echo 'specify variable FC via parmfile'
  exit 1
fi
export LD=$FC  # assume linker=compiler

if [[ $PARALLELL = true ]] ; then
  echo "choose FC compatible with MPI (mpifort or so)"
  exit 1
fi

if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report variables after analysing system-- ----------------"
   echo "----------------------------------------------------------------------"
   LIST="MAKE AR CPP FC LD "
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:10}
      echo -e "$X: $Y" 
   done
   echo "-------------------------------------------------------done-----------"
fi

################################################################################
##     process flags
##      -- attach parallel flag
##      -- prepend -D to all entries in CPPFLAGS
################################################################################
#      CPPVAR_FEAST
#      CPPVAR_JADAMILU
#      CPPVAR_SLEPC  (#INCLUDE <FINCLUDE/SLEPCEPSDEF.H>)

if [[ $PARALLEL = true && \
      -z $(echo $(CPPFLAGS) | grep "CPPVARIABLE_PARALLEL") ]] ; then
  CPPFLAGS='${CPPFLAGS} -DCPPVARIABLE_PARALLEL'
  echo "Warning: cpp flag "-DCPPVARIABLE_PARALLEL" required in parallel mode."
  echo "         Flag has been added to CPPFLAGS"
fi

#-------------------------------------------------------------------------------
#--------------report FLAGS----------------------------------------------------
#-------------------------------------------------------------------------------
if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report variables after setting flags----------------------"
   echo "----------------------------------------------------------------------"
   LIST="CPPFLAGS FFLAGS LDFLAGS"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:10}
      echo -e "$X: $Y" 
   done
   echo "-------------------------------------------------------done-----------"
fi

################################################################################
##     search Libraries: BLAS LAPACK FFTW LIBXC MPI
################################################################################
# there is a filesystem hierarchy Standard (FHS)

# The Filesystem Hierarchy Standard describes the filesystem conventions
# of a Linux system. In this standard, folders /lib, /usr/lib and
# /usr/local/lib are the default folders to store shared libraries. The
# /lib folder has libraries used during the boot time of the system but
# also used by programs in the /bin folder. Similarly, the/usr/lib
# folder has libraries used by programs in the /usr/bin folder. Finally,
# /usr/local/lib folder has libraries used by programs in /usr/local/bin
# folder.

# extension "so" (shared object) denotes a dynamic library
# extension "a" (archive) denotes a static library
# extension "dylib" (dynamic library) denotes a dynamic library on macOS
# https://stackoverflow.com/questions/2339679/what-are-the-differences-between-so-and-dylib-on-macos

if [[ $VERBOSE = true ]] ; then 
   echo  "searching libraries ............................................."
fi

#-------------------------------------------------------------------------------
#--- find libraries already on the libstring                                  --
#-------------------------------------------------------------------------------
DIR=
for X in $LIBS ; do
  if [[ ${X:0:2} = -L ]] ; then 
    DIR=${X#-L} 
  elif [[ ${X:0:2} = -l ]] ; then
    NAME=${X#-l}
    LIB=($(ls ${DIR}/lib${NAME}.*))
    if [[ -n $LIB ]] ; then
      case $NAME in
         blas)   BLASLIB=$LIB ;;
         lapack) LAPACKLIB=$LIB ;;
         fftw3)  FFTW3LIB=$LIB ;;
         xcf03)  XCLIB=$LIB ;;
         mpi)    MPILIB=$LIB ;;
      esac
    fi
  fi
done
if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------ libraries identified by parameter file-------------------"
   echo "----------------------------------------------------------------------"
   LIST="BLASLIB LAPACKLIB FFTW3LIB XCLIB MPILIB "
   for X in ${LIST}; do
      eval "Y=\${$X}"
      if [[ -z $Y ]] ; then continue ; fi
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:20}
      echo -e "$X: $Y" 
   done
   echo "-------------------------------------------------------done-----------"
fi

#-------------------------------------------------------------------------------
#--- search for missing libraries                                           ----
#-------------------------------------------------------------------------------
LIST=""
for X in BLASLIB LAPACKLIB FFTW3LIB XCLIB MPILIB ; do
   eval "Y=\${$X}"
  if [[ -z ${Y} ]] ; then LIST="$LIST $X" ; fi
done

# for macOS, remove blas and lapack. which are included via accelerate
if [[ $OS = Darwin ]] ; then 
  LIST=${LIST[@]/BLASLIB}
  LIST=${LIST[@]/LAPACKLIB}
fi

if [[ ${VERBOSE} = true ]] ; then
  echo "searching on the system for $LIST"
fi


# search through the standard directories with standard extensions
for X in $LIST ; do
   eval "Y=\${$X}"
   PATTERN=$(echo "$X" | tr '[:upper:]' '[:lower:]')
   PATTERN=lib${PATTERN%lib}
#   echo $PATTERN
   LOCS="/opt/homebrew /lib /usr/lib /usr/include /usr/local /opt"
   for LOC in $LOCS ; do
     if [[ ! -d $LOC ]] ; then continue; fi
     for EXT in dylib so a  ; do
       if [[ -n $Y ]] ; then continue ; fi # skip if already successful
       Z=($(find $LOC -name "$PATTERN*.$EXT" -print))
       if [[ ! -z $Z ]] ; then 
#        -- file encountered add to LIBS
         eval "$X=${Z[0]}"
         eval "Y=\${$X}"
         DIR=${Y%/*}
         NAME=${Y##*/lib}
         NAME=${NAME%%.*}
         LIBS="$LIBS -L$DIR -l$NAME"
#         echo X=$X $Y -L$DIR -l$NAME
         break 2   # break the loop over $LOCS
       fi 
     done
  done
done

#-------------------------------------------------------------------------------
#--- remove extra blanks from library string------------------------------------
#-------------------------------------------------------------------------------
LIBS=$(echo ${LIBS} | tr -s '[:blank:]')  # strip extra spaces

if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------libraries found on the system ----------------------------"
   echo "----------------------------------------------------------------------"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:20}
      echo -e "$X: $Y" 
   done
   echo ""
   echo "library string LIBS=$LIBS"
   echo "-------------------------------------------------------done-----------"
fi


################################################################################
##     Include and module files                                               ##
################################################################################
#-------------------------------------------------------------------------------
#-- identify files on the include string
#-------------------------------------------------------------------------------
LOCS="/lib /usr/lib /usr/local/lib /opt/homebrew opt"

#_________search for the includ file required for fftw3_________________________
if [[ -z $(echo $INCLUDES | grep 'fftw3.f03') ]] ; then
echo here
  for LOC in $LOCS ; do
   if [[ ! -d $LOC ]] ; then continue; fi
    X=($(find $LOC -name "fftw3.f03" -print))
    if [[ -n $X ]] ; then 
      X=${X[0]} 
      INCLUDES="$INCLUDES ${X}"
      break
    fi 
  done
fi

#_________search for the includ file required for libxc_________________________
if [[ -z $(echo $INCLUDES | grep 'xc_f03_lib_m.mod') ]] ; then
  for LOC in $LOCS ; do
   if [[ ! -d $LOC ]] ; then continue; fi
    X=($(find $LOC -name "xc_f03_lib_m.mod" -print))
    if [[ -n $X ]] ; then 
      X=${X[0]} 
      INCLUDES="$INCLUDES ${X}"
      break
    fi 
  done
fi

#_________search for the includ file required for mpi___________________________
if [[ -z $(echo $INCLUDES | grep 'mpi_f08.mod') ]] ; then
  for LOC in $LOCS ; do
   if [[ ! -d $LOC ]] ; then continue; fi
    X=($(find $LOC -name 'mpi_f08.mod' -print))
    if [[ -n $X ]] ; then 
      X=${X[0]} 
      INCLUDES="$INCLUDES ${X}"
      break
    fi 
  done
fi

#-------------------------------------------------------------------------------
#--- remove extra blanks from include string------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#--- report include files                                                     --
#-------------------------------------------------------------------------------
if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------include files---------------------------------------------"
   echo "----------------------------------------------------------------------"
   for X in ${INCLUDES}; do
      echo -e "$X" 
   done
   echo "include string INCLUDES=$INCLUDES"
   echo "-------------------------------------------------------done-----------"
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
##  report settings after fix
################################################################################
if [[ $VERBOSE = true ]] ; then
   echo "----------------------------------------------------------------------"
   echo "------------report variables buildfix---------------------------------"
   echo "----------------------------------------------------------------------"
  LIST="SUFFIX PARALLEL MAKE AR CPP FC LD \
        CPPFLAGS FFLAGS LDFLAGS \
        INCLUDES LIBS \
        BASEDIR BUILDDIR BINDIR"
   for X in ${LIST}; do
      eval "Y=\${$X}"
      Y=$(echo $Y | tr -s '[:blank:]')  # strip extra spaces
      X="${X}.............................."
      X=${X:0:20}
      echo -e "$X: $Y" 
   done
   echo "-------------------------------------------------------done-----------"
fi
