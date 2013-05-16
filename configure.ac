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
if test "$PARMFILE" = "notset"; then
 AC_MSG_ERROR(parmfile must be specified: use option --with-parmfile=parmfile)
else
 AC_CHECK_FILE($PARMFILE,TCHK="T",TCHK="F")
 if test "$TCHK" = "F"; then
   AC_MSG_ERROR(parmfile $PARMFILE not found)
 fi
fi
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
##                  check for fortran compilers                           ##
##            currently not used                                          ##
############################################################################
# g95  compiler (Open source)
AC_CHECK_PROG(EXIST_COMPILER_G95,g95,yes,no)
# xlf90 compiler (IBM)
AC_CHECK_PROG(EXIST_COMPILER_XLF90,xlf90,yes,no)
# ifort/ifc (Intel)
AC_CHECK_PROG(EXIST_COMPILER_IFORT,ifort,yes,no)
AC_CHECK_PROG(EXIST_COMPILER_IFC,ifc,yes,no)
# f90 (Absoft)
AC_CHECK_PROG(EXIST_COMPILER_F90,f90,yes,no)
# fort (F90)
AC_CHECK_PROG(EXIST_COMPILER_FORT,fort,yes,no)
#
############################################################################
##                                                                        ##
##                  check for libraries                                   ##
##               currently not used                                       ##
############################################################################
#
############################################################################
##                                                                        ##
##                  resolve parameter file                                ##
##                                                                        ##
############################################################################
#
AC_MSG_RESULT(resolve $PARMFILE)
if test "$PARMFILE" != notset; then
  ARCH=`grep '^ARCH' $PARMFILE | cut -d'"' -f2`
  TPARALLEL=`grep '^TPARALLEL' $PARMFILE | cut -d'"' -f2`
  TUPPERCASEMOD=`grep '^TUPPERCASEMOD' $PARMFILE | cut -d'"' -f2`
  COMPILER_SCALAR=`grep '^COMPILER_SCALAR' $PARMFILE | cut -d'"' -f2`
  COMPILER_PARALLEL=`grep '^COMPILER_PARALLEL' $PARMFILE | cut -d'"' -f2`
  FCFLAGS_NONE=`grep '^FCFLAGS_NONE' $PARMFILE | cut -d'"' -f2`
  FCFLAGS_OPT=`grep '^FCFLAGS_OPT' $PARMFILE | cut -d'"' -f2`
  FCFLAGS_DBG=`grep '^FCFLAGS_DBG' $PARMFILE | cut -d'"' -f2`
  FCFLAGS_PROF=`grep '^FCFLAGS_PROF' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_NONE=`grep '^LDFLAGS_NONE' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_OPT=`grep '^LDFLAGS_OPT' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_DBG=`grep '^LDFLAGS_DBG' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_PROF=`grep '^LDFLAGS_PROF' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_SCALAR=`grep '^LDFLAGS_SCALAR' $PARMFILE | cut -d'"' -f2`
  LDFLAGS_PARALLEL=`grep '^LDFLAGS_PARALLEL' $PARMFILE | cut -d'"' -f2`
  LIBS_SCALAR=`grep '^LIBS_SCALAR=' $PARMFILE | cut -d'"' -f2`
  LIBS_PARALLEL=`grep '^LIBS_PARALLEL' $PARMFILE | cut -d'"' -f2`
  CPPFLAGS=`grep '^CPPFLAGS' $PARMFILE | cut -d'"' -f2`
  FEXT=`grep '^FEXT' $PARMFILE | cut -d'"' -f2`
  BLASDIR=`grep '^BLASDIR' $PARMFILE | cut -d'"' -f2`
  LAPACKDIR=`grep '^LAPACKDIR' $PARMFILE | cut -d'"' -f2`
  FFTDIR=`grep '^FFTDIR' $PARMFILE | cut -d'"' -f2`
  MPIDIR=`grep '^MPIDIR' $PARMFILE | cut -d'"' -f2`
  FFT_HEADER=`grep '^FFT_HEADER' $PARMFILE | cut -d'"' -f2`
  MPI_HEADER=`grep '^MPI_HEADER' $PARMFILE | cut -d'"' -f2`
  SPECIAL=`grep '^SPECIAL' $PARMFILE | cut -d'"' -f2`
  #  
  # compose final variables
  #
  PAWDIR=`pwd`  # This is for absolute paths  
#  OBJDIR="$OBJDIR/$ARCH/Objects"
  PAWX="$PAWDIR/bin/$ARCH"
  COMPILE_SCALAR="$COMPILER_SCALAR"
  COMPILE_PARALLEL="$COMPILER_PARALLEL"
  LINK_SCALAR="$COMPILER_SCALAR $LDFLAGS_SCALAR"
  LINK_PARALLEL="$COMPILER_PARALLEL $LDFLAGS_PARALLEL"
  CPPFLAGS="$CPPFLAGS"
  FEXT="$FEXT"
  ARCH="$ARCH"
  if test "$TUPPERCASEMOD" != "T"; then 
    TUPPERCASEMOD="F"
  fi
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
echo '${MPI_HEADER}'
AC_DEFINE([XX],${MPI_HEADER})
echo ${XX}
echo '============================='
if test ["$TPARALLEL" = "T"]; then
  #________mpi include file
  if test [ -n "$MPIDIR" ] ; then
    AC_MSG_RESULT(check MPI directory)
    AC_CHECK_FILES($MPIDIR,ERROR=F,ERROR=T)
    if test ["$ERROR" = "T"]; then
      AC_MSG_RESULT(name of the mpi directory: $MPIDIR)
      AC_MSG_ERROR(mpi directory not found)
    fi
##    AC_MSG_RESULT(check MPI include file)
##    AC_CHECK_FILES($MPI_HEADER,ERROR=F,ERROR=T)
##    if test ["$ERROR" = "T"]; then
##      AC_MSG_RESULT(name of the mpi include file: $MPIDIR)
##      AC_MSG_ERROR(mpi include file not found)
##    fi
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
##  AC_MSG_RESULT(check FFT include file)
##  AC_CHECK_FILES($FFT_HEADER,ERROR=F,ERROR=T)
##  if test ["$ERROR" = "T"]; then
##    AC_MSG_RESULT(name of the fftw include file: $FFT_HEADER)
##    AC_MSG_ERROR(fftw include file not found)
##  fi
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
AC_SUBST(COMPILE_SCALAR)
AC_SUBST(COMPILE_PARALLEL)
AC_SUBST(LINK_SCALAR)
AC_SUBST(LINK_PARALLEL)
AC_SUBST(LIBS_SCALAR)
AC_SUBST(LIBS_PARALLEL)
AC_SUBST(FEXT)
AC_SUBST(CPPFLAGS)
AC_SUBST(PAWX)     
AC_SUBST(TPARALLEL)
#
AC_SUBST(BLASDIR)     
AC_SUBST(LAPACKDIR)     
AC_SUBST(FFTDIR)     
AC_SUBST(MPIDIR)     
AC_SUBST(FFT_HEADER)     
AC_SUBST(MPI_HEADER)     
#
#########################################################################
##                                                                     ##
##  create Makefiles and f90pp                                         ##
##                                                                     ##
#########################################################################
#
# copy current parameter file to "parms.in_use"
#
AC_MSG_RESULT(copy $PARMFILE to parms.in_use)
cp $PARMFILE ${PAWDIR}/parms.in_use
#
# creating subdirectories
#
AC_MSG_RESULT(creating subdirectories and copying shell scripts)
mkdir -p  ${PAWDIR}/bin/${ARCH} 
mkdir -p  ${PAWDIR}/doc
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
# substituting to accept lowercase module files
# default are uppercase module file names with lowercase extension
#
if test "$TUPPERCASEMOD" = "F"; then
  AC_MSG_RESULT(modify makefile_bare for lowercase module files)
#  sed -e 's%\([[:space]].*.mod\)%\L\1%'\
# pattern matching: [^\/] do not match a slash; 
#                   (/string/) save string for later use as /1
#                   \L make expression lowercase
# sed command       s%stringa%stringb%g  replace stringa by stringb 
#
#  sed -e 's%\([[^\/]]*.mod\)%\L\1%'\
#      ${PAWDIR}/Makefile_bare > ${PAWDIR}/Makefile_bare_tmp
  ./lowercasemodnames.sh  < ${PAWDIR}/Makefile_bare > ${PAWDIR}/Makefile_bare_tmp
  cp  ${PAWDIR}/Makefile_bare_tmp  ${PAWDIR}/Makefile_bare
  rm  -f ${PAWDIR}/Makefile_bare_tmp
fi
#
# substituting to get Makefile.none
#
AC_MSG_RESULT(creating Makefile in none)
sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_NONE %g"\
    -e "s%\(LINK=.*\)%\1 $LDFLAGS_NONE %g"\
    -e "s%\(PAWX.*=.*\)%\1/paw.x%g"\
    -e 's%varPARALLEL%F%g'\
    -e 's%varOBJTYPE%none%g'\
    ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/none/Makefile
#
# substituting to get Makefile_parallel in none
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile in none_parallel)
  sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_NONE %g"\
      -e "s%\(LINK=.*\)%\1 $LDFLAGS_NONE %g"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw.x%g"\
      -e 's%varPARALLEL%T%g'\
      -e 's%varOBJTYPE%none_parallel%g'\
      ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/none_parallel/Makefile
fi
#
# substituting to get Makefile.fast
#
AC_MSG_RESULT(creating Makefile in fast)
sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_OPT %g"\
    -e "s%\(LINK=.*\)%\1 $LDFLAGS_OPT %g"\
    -e "s%\(PAWX.*=.*\)%\1/paw_fast.x%g"\
    -e 's%varOBJTYPE%fast%g'\
    -e 's%varPARALLEL%F%g'\
    ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/fast/Makefile
#
# substituting to get Makefile.fast_parallel
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in fast)
  sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_OPT %g"\
      -e "s%\(LINK=.*\)%\1 $LDFLAGS_OPT %g"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_fast.x%g"\
      -e 's%varOBJTYPE%fast_parallel%g'\
      -e 's%varPARALLEL%T%g'\
      ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/fast_parallel/Makefile
fi
#
# substituting to get Makefile in dbg
#
AC_MSG_RESULT(creating Makefile in dbg)
sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_DBG %g"\
    -e "s%\(LINK=.*\)%\1 $LDFLAGS_DBG %g"\
    -e "s%\(PAWX.*=.*\)%\1/paw_dbg.x%g"\
    -e 's%varOBJTYPE%dbg%g'\
    -e 's%varPARALLEL%F%g'\
   ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/dbg/Makefile
#
# substituting to get Makefile_parallel in dbg
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in dbg)
  sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_DBG %g"\
      -e "s%\(LINK=.*\)%\1 $LDFLAGS_DBG %g"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_dbg.x%g"\
      -e 's%varOBJTYPE%dbg_parallel%g'\
      -e 's%varPARALLEL%T%g'\ 
     ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/dbg_parallel/Makefile
fi
#
# substituting to get Makefile in prof
#
AC_MSG_RESULT(creating Makefile in prof)
sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_PROF %g"\
    -e "s%\(LINK=.*\)%\1 $LDFLAGS_PROF %g"\
    -e "s%\(PAWX.*=.*\)%\1/paw_prof.x%g"\
    -e 's%varOBJTYPE%prof%g'\
    -e 's%varPARALLEL%F%g'\
     ${PAWDIR}/Makefile_bare > ${PAWDIR}/bin/${ARCH}/Objects/prof/Makefile
#
# substituting to get Makefile_parallel in prof
#
if test "$TPARALLEL" = "T"; then
  AC_MSG_RESULT(creating Makefile_parallel in prof)
  sed -e "s%\(COMPILE=.*\)%\1 $FCFLAGS_PROF %g"\
      -e "s%\(LINK=.*\)%\1 $LDFLAGS_PROF %g"\
      -e "s%\(PAWX.*=.*\)%\1/ppaw_prof.x%g"\
      -e 's%varOBJTYPE%prof_parallel%g'\
      -e 's%varPARALLEL%T%g'\
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
AC_MSG_RESULT(---------------------------------------------------------------)
AC_MSG_RESULT(-----------------------SUMMARY---------------------------------)
AC_MSG_RESULT(---------------------------------------------------------------)
AC_MSG_RESULT(directory of distribution      : $PAWDIR)
AC_MSG_RESULT(directory with binaries        : $PAWDIR/bin/$ARCH)
AC_MSG_RESULT(architecture                   : $ARCH)
AC_MSG_RESULT(preprocessor variables         : $CPPFLAGS)
AC_MSG_RESULT(architecture name              : $ARCH)
AC_MSG_RESULT(parallel environment           : $TPARALLEL)
AC_MSG_RESULT(compile command (scalar)       : $COMPILE_SCALAR)
AC_MSG_RESULT(F90 file extension             : $FEXT)
AC_MSG_RESULT(compile flags (none)           : $FCFLAGS_NONE)
AC_MSG_RESULT(compile flags (fast)           : $FCFLAGS_OPT)
AC_MSG_RESULT(compile flags (dbg)            : $FCFLAGS_DBG)
AC_MSG_RESULT(compile flags (prof)           : $FCFLAGS_PROF)
AC_MSG_RESULT(special link flags (none)      : $LDFLAGS_NONE)
AC_MSG_RESULT(special link flags (opt)       : $LDFLAGS_OPT)
AC_MSG_RESULT(special link flags (dbg)       : $LDFLAGS_DBG)
AC_MSG_RESULT(special link flags (prof)      : $LDFLAGS_PROF)
AC_MSG_RESULT(link command w.flags (scalar)  : $LINK_SCALAR)
AC_MSG_RESULT(external libraries (scalar)    : $LIBS_SCALAR)
AC_MSG_RESULT(uppercase module names?        : $TUPPERCASEMOD)
#AC_MSG_RESULT(---------------------------------------------------------------)
#AC_MSG_RESULT(--------------external libraries-------------------------------)
#AC_MSG_RESULT(---------------------------------------------------------------)
AC_MSG_RESULT(blas library                   : $BLASDIR)
AC_MSG_RESULT(lapack library                 : $LAPACKDIR)
AC_MSG_RESULT(libs for Fourier transforms    : $FFTDIR)
if test "$TPARALLEL" = "T"; then
#  AC_MSG_RESULT(---------------------------------------------------------------)
#  AC_MSG_RESULT(--------------specifics for parallel environment---------------)
#  AC_MSG_RESULT(---------------------------------------------------------------)
  AC_MSG_RESULT(parallel envirnment considered : yes)
  AC_MSG_RESULT(compile command (parallel)     : $COMPILE_PARALLEL)
  AC_MSG_RESULT(link command w. flags(parallel): $LINK_PARALLEL)
  AC_MSG_RESULT(external libraries (parallel)  : $LIBS_PARALLEL)
  AC_MSG_RESULT(mpi library                    : $MPIDIR)
else
  AC_MSG_RESULT(parallel envirnment considered : no)
fi
#AC_MSG_RESULT(---------------------------------------------------------------)
#AC_MSG_RESULT(--------------environment related------------------------------)
#AC_MSG_RESULT(---------------------------------------------------------------)
AC_MSG_RESULT(MAKE command                   : $MAKE)
AC_MSG_RESULT(CPP command                    : $CPP)
AC_MSG_RESULT(---------------------------------------------------------------)
AC_MSG_RESULT(----done!---configuration completed successfully!--------------)
AC_MSG_RESULT(---------------------------------------------------------------)
