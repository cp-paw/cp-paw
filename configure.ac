dnl Process this file with autoconf to produce a configure script.
#
#
# This line checks if the source-directory is in the actual path
# This prevents users to call the configure script from somewhere else
## MIND !!! adapt next line if changing the name of the source directory
#
AC_INIT(src)
#
#######################################################################
#######################################################################
##                                                                   ##
## begin of user adaptable variables                                 ##
##                                                                   ##
#######################################################################
#######################################################################
#
# directory structure is set to
#   $PAWDIR/src               contains source files for simulation code
#   $PAWDIR/F90PP/            contains files for f90pp preprocessor
#   $PAWDIR/Tools/            contains subdirectories with tools
#   $PAWDIR/bin/$ARCH         holds all executables
#   $PAWDIR/bin/$ARCH/none    Objects compiled simple
#   $PAWDIR/bin/$ARCH/prof    Objects compiled for profiling
#   $PAWDIR/bin/$ARCH/dbg     Objects compiled with debugger
#   $PAWDIR/bin/$ARCH/fast    Objects compiled with optimizer
#
# PAWDIR will be set automatically after this block!
#
#_________________way to compile the PAW-tools
TOOLCOMP="dbg"         
#
#######################################################################
#######################################################################
##                                                                   ##
## end of user adaptable variables                                   ##
##                                                                   ##
#######################################################################
#######################################################################
#
# subdirectory tree of PAW distribution
#
# the PAWDIR is determined automatically, the rest follows 
#
PAWDIR=`pwd`
#
# process command line arguments
#
# syntax: AC_ARG_WITH(package,help string,action if given, action if not given)
# can be invoked by calling ./configure --with-package= or
# ./configure --without-package
#
AC_ARG_WITH(make,--with-make=LOCATION\
      the location of the GNU make utility (if not in path),
            MAKEPATH=$withval,MAKEPATH=notset)
AC_ARG_WITH(cpp,--with-cpp=LOCATION\
       the location of the GNU C Preprocessor (if not in path),
            CPPPATH=$withval,CPPPATH=notset)
AC_ARG_WITH(parmfile,--with-parmfile=PARMFILE\
  specify parameterfile - no automatic configuration,
            PARMFILE=$withval,PARMFILE=notset)
#
# Check for GNU make utility
#
# set by command line argument?
if test "$MAKEPATH" = "notset"; then
  # syntax: AC_CHECK_PROGS(variable, programs, value if found, value if not found)
  AC_MSG_RESULT(check for make)
  AC_CHECK_PROGS(MAKE, gmake make,notset)
else
  AC_MSG_RESULT(----------------------------------------------)
  AC_MSG_RESULT(Maketool set by command line agrument to)
  AC_MSG_RESULT($MAKEPATH)
  AC_MSG_RESULT(----------------------------------------------)
  AC_CHECK_FILE($MAKEPATH,MAKE="$MAKEPATH",MAKE="notset")
fi

if test "$MAKE" = "notset"; then
   AC_MSG_ERROR(GNU make utility was not found)
fi

#
# Check for the GNU C-Compatible Compiler Preprocessor
#
# set by command line argument?
if test "$CPPPATH" = "notset"; then
  AC_CHECK_PROG(CPP,cpp,cpp,notset)
else
  AC_MSG_RESULT(----------------------------------------------)
  AC_MSG_RESULT(C-Preprocessor set by command line argument to)
  AC_MSG_RESULT($CPPPATH)
  AC_MSG_RESULT(----------------------------------------------)
  AC_CHECK_FILE($CPPPATH,CPP="$CPPPATH",CPP="notset")
fi

if test "$CPP" = "notset"; then
   AC_MSG_ERROR(GNU C-Compatible Compiler Preprocessor is not found)
fi
CPP="$CPP $CPPOPT"

############################################################################
##                                                                        ##
##                  resolve parameter file                                ##
##                                                                        ##
############################################################################
#
# copy current parameter file to "parms.in_use"
#
AC_MSG_RESULT(copy $PARMFILE to parms.in_use)
cp $PARMFILE ${PAWDIR}/parms.in_use
#
AC_MSG_RESULT(resolve $PARMFILE)
if test "$PARMFILE" != notset; then
  ARCH=`grep '^ARCH' $PARMFILE | cut -d'"' -f2`
  TPARALLEL=`grep '^TPARALLEL' $PARMFILE | cut -d'"' -f2`
  COMPILER=`grep '^COMPILER' $PARMFILE | cut -d'"' -f2`
  FLAGS=`grep '^FLAGS' $PARMFILE | cut -d'"' -f2`
  OPTIMIZE=`grep '^OPTIMIZE' $PARMFILE | cut -d'"' -f2`
  DEBUG=`grep '^DEBUG' $PARMFILE | cut -d'"' -f2`
  PROF=`grep '^PROF' $PARMFILE | cut -d'"' -f2`
  LFLAGS=`grep '^LFLAGS' $PARMFILE | cut -d'"' -f2`
  LLINK=`grep '^LLINK' $PARMFILE | cut -d'"' -f2`
  PLLINK=`grep '^PLLINK' $PARMFILE | cut -d'"' -f2`
  LIBS=`grep '^LIBS=' $PARMFILE | cut -d'"' -f2`
  PLIBS=`grep '^PLIBS' $PARMFILE | cut -d'"' -f2`
  CPPVAR=`grep '^CPPVAR' $PARMFILE | cut -d'"' -f2`
  FEXT=`grep '^FEXT' $PARMFILE | cut -d'"' -f2`
  BLASDIR=`grep '^BLASDIR' $PARMFILE | cut -d'"' -f2`
  LAPACKDIR=`grep '^LAPACKDIR' $PARMFILE | cut -d'"' -f2`
  FFTDIR=`grep '^FFTDIR' $PARMFILE | cut -d'"' -f2`
  MPIDIR=`grep '^MPIDIR' $PARMFILE | cut -d'"' -f2`
  INCLUDE_FFTW=`grep '^INCLUDE_FFTW' $PARMFILE | cut -d'"' -f2`
  INCLUDE_MPI=`grep '^INCLUDE_MPI' $PARMFILE | cut -d'"' -f2`
  SPECIAL=`grep '^SPECIAL' $PARMFILE | cut -d'"' -f2`
  #  
  # compose final variables
  #
  COMPILE="$COMPILER $FLAGS"
  LINK="$COMPILER $LFLAGS"
  LLINK="$LLINK"
  LLINKPARALLEL="$PLLINK"
  LLIBS="$LIBS"
  LLIBSPARALLEL="$PLIBS"
  LLIBSPARALLEL="$LLIBSPARALLEL $USERLIBS"
  CPPVAR="$CPPVAR"
  FEXT="$FEXT"
  ARCH="$ARCH"
  OBJDIR="$OBJDIR/$ARCH/Objects"
  PAWX="$PAWDIR/bin/$ARCH"
  #
  # if no parameter file is given, the script will search for the
  # values needed.
  #
else
   AC_MSG_ERROR(No Parameter file. Specify --with-parmfile=parameterfile)
fi     # end of block "without parmfile" #--------------------------------------

#
# Definition of variables to be substituted in the *.in files
#   There the synax is @VARIABLE@
#
#ATTENTION CHANGE!!!!
CPP="$CPP -traditional"
#
#########################################################################
##                                                                     ##
##  check if files are present                                         ##
##                                                                     ##
#########################################################################
if test ["$TPARALLEL" = "T"]; then
  #________mpi include file
  if test [ -n "$MPIDIR" ] ; then
    AC_MSG_RESULT(check MPI directory)
    AC_CHECK_FILES($MPIDIR,ERROR=F,ERROR=T)
    if test ["$ERROR" = "T"]; then
      AC_MSG_RESULT(name of the mpi directory: $MPIDIR)
      AC_MSG_ERROR(mpi directory not found)
    fi
    AC_MSG_RESULT(check MPI include file)
    AC_CHECK_FILES($INCLUDE_MPI,ERROR=F,ERROR=T)
    if test ["$ERROR" = "T"]; then
      AC_MSG_RESULT(name of the mpi include file: $MPIDIR)
      AC_MSG_ERROR(mpi include file not found)
    fi
  else
    # mpidir is empty. error or ok?
    if test "$SPECIAL"="IBM" ; then
      AC_MSG_RESULT(POE environment: mpi directory not needed)
    else
      AC_MSG_ERROR(mpi directory not specified)
    fi
  fi
fi
#
#  fft directory
if test -n "$FFTDIR" ; then
  AC_MSG_RESULT(check FFT directory)
  AC_CHECK_FILES($FFTDIR,ERROR=F,ERROR=T)
  if test ["$ERROR" = "T"]; then
    AC_MSG_RESULT(name of the FFT directory: $FFTDIR)
    AC_MSG_ERROR(FFT directory not found)
  fi
  AC_MSG_RESULT(check FFTW include file)
  AC_CHECK_FILES($INCLUDE_FFTW,ERROR=F,ERROR=T)
  if test ["$ERROR" = "T"]; then
    AC_MSG_RESULT(name of the fftw include file: $INCLUDE_FFTW)
    AC_MSG_ERROR(fftw include file not found)
  fi
fi
#
#########################################################################
##                                                                     ##
##  define substitutions for AC_OUTPUT                                 ##
##                                                                     ##
#########################################################################
AC_SUBST(MAKE)
AC_SUBST(CPP)
AC_SUBST(PAWDIR)
AC_SUBST(ARCH)
AC_SUBST(SPECIAL)
AC_SUBST(COMPILE)
AC_SUBST(LINK)
AC_SUBST(LLINK)
AC_SUBST(LLINKPARALLEL)
AC_SUBST(LLIBSPARALLEL)
AC_SUBST(LLIBS)
AC_SUBST(FEXT)
AC_SUBST(TOOLCOMP)
AC_SUBST(CPPVAR)
AC_SUBST(PAWX)     
#
AC_SUBST(BLASDIR)     
AC_SUBST(LAPACKDIR)     
AC_SUBST(FFTDIR)     
AC_SUBST(MPIDIR)     
AC_SUBST(INCLUDE_FFTW)     
AC_SUBST(INCLUDE_MPI)     
#
#########################################################################
##                                                                     ##
##  create Makefiles and f90pp                                         ##
##                                                                     ##
#########################################################################
#
# creating subdirectories
#
AC_MSG_RESULT(creating subdirectories and copying shell scripts)
mkdir -p  ${PAWDIR}/bin/${ARCH} 
mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/none
mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/dbg
mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/prof
mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/fast
if test "$TPARALLEL" = "T"; then
  mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/none_parallel
  mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/dbg_parallel
  mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/prof_parallel
  mkdir -p  ${PAWDIR}/bin/${ARCH}/Objects/fast_parallel
fi
#
#  copy script files from
#
#AC_MSG_RESULT(copy scripts into $PAWDIR/bin/$ARCH/)
#cp $PAWDIR/src/Tools/Scripts/* $PAWDIR/bin/$ARCH/
#
# create f90pp from f90pp.in and make it executable
#
AC_CONFIG_FILES($PAWDIR/bin/$ARCH/f90pp:$PAWDIR/src/F90PP/f90pp.in,[chmod +x $PAWDIR/bin/$ARCH/f90pp])
#
# create Makefile from Makefile_targets.in 
#
AC_CONFIG_FILES(Makefile:Makefile_targets.in)
#
# create Makefile_bare from Makefile.in 
#
AC_CONFIG_FILES(Makefile_bare:Makefile_bare.in)
#
#  NOW OUTPUT THE FILES CREATED BY ACCONFIG_FILES
#
AC_OUTPUT()
chmod +x $PAWDIR/bin/$ARCH/f90pp
#
# substituting to get Makefile.none
#
AC_MSG_RESULT(creating Makefile in none)
sed -e "s%\(PAWX.*=.*\)%\1/paw.x%"\
    -e 's%varPARALLEL%F%'\
    -e 's%varOBJTYPE%none%'\
    ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/none/Makefile
#
# substituting to get Makefile_parallel in none
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile in none_parallel)
  sed -e "s%\(PAWX.*=.*\)%\1/ppaw.x%"\
      -e 's%varPARALLEL%T%'\
      -e 's%varOBJTYPE%none_parallel%'\
      ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/none_parallel/Makefile
fi
#
# substituting to get Makefile.fast
#
AC_MSG_RESULT(creating Makefile in fast)
sed -e "s%\(COMPILE=.*\)%\1 $OPTIMIZE %"\
    -e "s%\(PAWX.*=.*\)%\1/paw_fast.x%"\
    -e 's%varOBJTYPE%fast%'\
    -e 's%varPARALLEL%F%'\
    ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/fast/Makefile
#
# substituting to get Makefile.fast_parallel
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in fast)
  sed -e "s%\(COMPILE=.*\)%\1 $OPTIMIZE %"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_fast.x%"\
      -e 's%varOBJTYPE%fast_parallel%'\
      -e 's%varPARALLEL%T%'\
      ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/fast_parallel/Makefile
fi
#
# substituting to get Makefile in dbg
#
AC_MSG_RESULT(creating Makefile in dbg)
sed -e "s%\(COMPILE=.*\)%\1 $DEBUG %"\
    -e "s%\(PAWX.*=.*\)%\1/paw_dbg.x%"\
    -e 's%varOBJTYPE%dbg%'\
    -e 's%varPARALLEL%F%'\
   ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/dbg/Makefile
#
# substituting to get Makefile_parallel in dbg
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in dbg)
  sed -e "s%\(COMPILE=.*\)%\1 $DEBUG %"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_dbg.x%"\
      -e 's%varOBJTYPE%dbg_parallel%'\
      -e 's%varPARALLEL%T%'\ 
     ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/dbg_parallel/Makefile
fi
#
# substituting to get Makefile in prof
#
AC_MSG_RESULT(creating Makefile in prof)
sed -e "s%\(COMPILE=.*\)%\1 $PROF %"\
    -e "s%\(LINK=.*\)%\1 -pg %"\
    -e "s%\(PAWX.*=.*\)%\1/paw_prof.x%"\
    -e 's%varOBJTYPE%prof%'\
    -e 's%varPARALLEL%F%'\
     ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/prof/Makefile
#
# substituting to get Makefile_parallel in prof
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in prof)
  sed -e "s%\(COMPILE=.*\)%\1 $PROF %"\
      -e "s%\(LINK=.*\)%\1 -pg %"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_prof.x%"\
      -e 's%varOBJTYPE%prof_parallel%'\
      -e 's%varPARALLEL%T%'\
     ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/prof_parallel/Makefile
fi
#
# removal of config.cache 
#
#rm -f config.cache
#===== now output results
#########################################################################
##                                                                     ##
##                      print summary                                  ##
##                                                                     ##
#########################################################################
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(printing out of summary              )
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(architecture                   : $ARCH)
AC_MSG_RESULT(MAKE command                   : $MAKE)
AC_MSG_RESULT(CPP command                    : $CPP)
AC_MSG_RESULT(preprocessor variables         : $CPPVAR)
AC_MSG_RESULT(directory of distribution      : $PAWDIR)
AC_MSG_RESULT(architecture name              : $ARCH)
AC_MSG_RESULT(parallel environment           : $TPARALLEL)
AC_MSG_RESULT(F90 compile command            : $COMPILE)
AC_MSG_RESULT(F90 optimize options           : $OPTIMIZE)
AC_MSG_RESULT(F90 debug flags                : $DEBUG)
AC_MSG_RESULT(F90 profiling flags            : $PROF)
AC_MSG_RESULT(F90 linking flags              : $LFLAGS)
AC_MSG_RESULT(search directories (sequential): $LLINK)
AC_MSG_RESULT(search directories (parallel)  : $LLINKPARALLEL)
AC_MSG_RESULT(external libraries (sequential): $LLIBS)
AC_MSG_RESULT(external libraries (parallel)  : $LLIBSPARALLEL)
AC_MSG_RESULT(F90 file extension             : $FEXT)
AC_MSG_RESULT(blas library                   : $BLASDIR)
AC_MSG_RESULT(lapack library                 : $LAPACKDIR)
AC_MSG_RESULT(libs for Fourier transforms    : $FFTDIR)
AC_MSG_RESULT(mpi library                    : $MPIDIR)
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(configuration completed successfully!)
