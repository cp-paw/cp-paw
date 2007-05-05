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
# the default path for the binaries is
# $PAWDIR/$BINARIES/$*ARCH
# PAWDIR will be set automatically after this block!
# ________________names for the binaries
NONEBIN="paw.x"
FASTBIN="paw_fast.x"    
PROFBIN="paw_prof.x"    
DEBUGBIN="paw_dbg.x"  
PARALLELBIN="ppaw_fast.x"
#
# ________________directories for the architecture dependent binaries
#++++++absARCH="intel"   
#++++++fortARCH="alpha"   
xlf90ARCH="ibm"  
#------#
#------# Specification of the Fortran Compilers
#------# Fortran77 is needed for the internal numerical
#------# libraries
#------# ________________Fortran77/90 compilers 
#------AbsoftF90="f90"
#------AbsoftF77="f90"
#------CompaqF90="fort"
#------CompaqF77="fort"
#------IBMF90="xlf90"
#------IBMF77="xlf"
#------# ________________compiler flags for unoptimized compilation (Absoft,Compaq,IBM)
#------absFLAGS="-c -YEXT_NAMES=LCS -YEXT_SFX=_"
#------absf77FLAGS="$absFLAGS"
#------fortFLAGS="-c"
#------fort77FLAGS="$fortFLAGS"
#------xlf90FLAGS="-qlanglvl=90STD -c -qhalt=e -qextchk"
#------xlfFLAGS="-qlanglvl=77STD -c -qhalt=e -qextchk"
#------#
#------# ________________compiler flags for optimized compilation (Absoft,Compaq,IBM)
#------absOPTIMIZE90="-B100 -O"
#------fortOPTIMIZE90="-O3"
#------xlfOPTIMIZE90="-O3 -Q+clebsch -qfloat=rsqrt:hssngl:rsqrt -qmaxmem=-1 -qarch=auto"
#------absOPTIMIZE77="-B100 -O"
#------fortOPTIMIZE77="-O3"
#------xlfOPTIMIZE77="-O3 -Q+clebsch -qfloat=rsqrt:hssngl:rsqrt -qmaxmem=-1 -qarch=auto"
#------#
#------# ________________compiler flags for debug compilation (Absoft,Compaq,IBM)
#------absDEBUG90="-g -s"   
#------fortDEBUG90="-g3 -ladebug -C -assume gfullpath -check overflow -check underflow"
#------xlfDEBUG90="-g -C -qfullpath -qflttrap=zerodivide:invalid:enable -qsigtrap -qinitauto=FF"
#------absDEBUG77="-g -s"   
#------fortDEBUG77="-g3 -ladebug -C -assume gfullpath -check overflow -check underflow"
#------xlfDEBUG77="-g -C -qfullpath -qflttrap=zerodivide:invalid:enable -qsigtrap -qinitauto=FF"
#
# ________________compiler flags for linking the object files
#------absLFLAGS=""
#------fortLFLAGS=""
#------xlf90LFLAGS="-bmaxdata:0x70000000 -bloadmap:loadmap -qextchk"
#
# ________________source file extensions
#------absFEXT="f90"
#------fortFEXT="f90"
#------xlf90FEXT="f"
#
# ________________remove preprocessed sources after compilation?
#                 set to "rm" if yes and to "file" if no
NONEOPTRM="rm"
FASTOPTRM="rm"
PARALLELOPTRM="rm"
DEBUGOPTRM="file"
PROFOPTRM="file"
#
# the object files can be found in
# $PAWDIR/$BINARIES/$ARCH/$OBJECTS/$*DIR
# ________________directories for the objects
OBJECTS="Objects"
NONEDIR="none"  
FASTDIR="fast"
PROFDIR="prof"
DEBUGDIR="debug"       
PARALLELDIR="$FASTDIR" 
#_________________way to compile the PAW-tools
TOOLCOMP="debug"         
#
# the default path for the external libraries is
# $PAWDIR/$LIBRARIES/*ROOT
#
#________________ default directory for the libraries      
LIBRARIES="lib"  
#
#++++++#________________ default directory for FFTW-distribution
#++++++FFTWROOT="fftw"    # to be specified later! (fftw-$ARCH)
#++++++#
#++++++#________________ directory of MPICH-distribution
#++++++MPICHROOT="mpich"  # to be specified later! (mpich-$ARCH)
#++++++#
#++++++#________________ directory of ATLAS-distribution
#++++++ATLASROOT="ATLAS"  
#++++++#
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
SRCDIR="$PAWDIR/src" 
OBJDIR="$PAWDIR/bin"      # to be specified later (compiler-dependent)
BINDIR="$PAWDIR/bin"      # to be specified later (compiler-dependent) 
TOOLDIR="$PAWDIR/src/Tools"   
SCRIPTDIR="$PAWDIR/src/Tools/Scripts"   
#
# process command line arguments
#
# syntax: AC_ARG_WITH(package,help string,action if given, action if not given)
# can be invoked by calling ./configure --with-package= or
# ./configure --without-package
#
AC_ARG_WITH(libs,--without-libs\
            do NOT use external numerical libraries,
            USELIBS=$withval,USELIBS=yes)
AC_ARG_WITH(fftwdir,--with-fftwdir=DIR\
        root directory of FFTW distribution,
            FFTWDIR=$withval,FFTWDIR=notset)
AC_ARG_WITH(mpichdir,--with-mpichdir=DIR\
       root directory of MPICH distribution,
            MPICHDIR=$withval,MPICHDIR=notset)
#-----AC_ARG_WITH(atlasdir,--with-atlasdir=DIR\
#-----       root directory of ATLAS distribution,
#-----            ATLASDIR=$withval,ATLASDIR=notset)
#-----AC_ARG_WITH(atlasarch,--with-atlasarch=ARCH\
#-----     architecture for wich ATLAS blas was compiled (e.g. Linux_ATHLON),
#-----            ATLASARCH=$withval,ATLASARCH=notset)
AC_ARG_WITH(make,--with-make=LOCATION\
      the location of the GNU make utility (if not in path),
            MAKEPATH=$withval,MAKEPATH=notset)
AC_ARG_WITH(cpp,--with-cpp=LOCATION\
       the location of the GNU C Preprocessor (if not in path),
            CPPPATH=$withval,CPPPATH=notset)
AC_ARG_WITH(cppoptions,--with-cppoptions=OPTIONS\
 options for the GNU C Preprocessor,
            CPPOPT=$withval,CPPOPT="")
AC_ARG_WITH(cflags,--with-cflags=FLAGS\
       specify additional compiler flags,
            USERFLAGS=$withval)
AC_ARG_WITH(lflags,--with-lflags=LFLAGS\
      specify additional compiler flags for linking,
            USERLFLAGS=$withval)
AC_ARG_WITH(llibs,--with-llibs=LLIBS\
        specify additional libraries to be linked,
            USERLIBS=$withval)
AC_ARG_WITH(parmfile,--with-parmfile=PARMFILE\
  specify parameterfile - no automatic configuration,
            PARMFILE=$withval,PARMFILE=notset)

#
# Check for GNU make utility
#
# set by command line argument?
if test "$MAKEPATH" = "notset"; then
  # syntax: AC_CHECK_PROGS(variable, programs, value if found, value if not found)
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
if test "$PARMFILE" != notset; then
  ARCH=`grep '^ARCH' $PARMFILE | cut -d'"' -f2`
  COMPILER=`grep '^COMPILER' $PARMFILE | cut -d'"' -f2`
  F77COMPILER=`grep '^F77COMPILER' $PARMFILE | cut -d'"' -f2`
  FLAGS=`grep '^FLAGS' $PARMFILE | cut -d'"' -f2`
  F77FLAGS=`grep '^F77FLAGS' $PARMFILE | cut -d'"' -f2`
  OPTIMIZE=`grep '^OPTIMIZE' $PARMFILE | cut -d'"' -f2`
  F77OPTIMIZE=`grep '^F77OPTIMIZE' $PARMFILE | cut -d'"' -f2`
  DEBUG=`grep '^DEBUG' $PARMFILE | cut -d'"' -f2`
  F77DEBUG=`grep '^F77DEBUG' $PARMFILE | cut -d'"' -f2`
  PROF=`grep '^PROF' $PARMFILE | cut -d'"' -f2`
  F77PROF=`grep '^F77PROF' $PARMFILE | cut -d'"' -f2`
  LFLAGS=`grep '^LFLAGS' $PARMFILE | cut -d'"' -f2`
  LLINK=`grep '^LLINK' $PARMFILE | cut -d'"' -f2`
  PLLINK=`grep '^PLLINK' $PARMFILE | cut -d'"' -f2`
  LIBS=`grep '^LIBS=' $PARMFILE | cut -d'"' -f2`
  PLIBS=`grep '^PLIBS' $PARMFILE | cut -d'"' -f2`
  CPPVAR=`grep '^CPPVAR' $PARMFILE | cut -d'"' -f2`
  FEXT=`grep '^FEXT' $PARMFILE | cut -d'"' -f2`
  MPICHDIR=`grep '^MPICHDIR' $PARMFILE | cut -d'"' -f2`
  LAPACKDIR=`grep '^LAPACKDIR' $PARMFILE | cut -d'"' -f2`
  FFTWDIR=`grep '^FFTWDIR' $PARMFILE | cut -d'"' -f2`
  #  
  # compose final variables
  #
  COMPILE="$COMPILER $FLAGS"
  COMPILE77="$F77COMPILER $F77FLAGS"
  LINK="$COMPILER $LFLAGS"
  LLINK="$LLINK"
  LLINKPARALLEL="$PLLINK"
  LLIBS="$LIBS"
  LLIBSPARALLEL="$PLIBS"
  CPPVAR="$CPPVAR"
  FEXT="$FEXT"
  ARCH="$ARCH"
  OBJDIR="$OBJDIR/$ARCH/$OBJECTS"
  PAWX="$BINDIR/$ARCH"
  #
  # if no parameter file is given, the script will search for the
  # values needed.
  #
else
   AC_MSG_ERROR(No Parameter file. Specify --with-parmfile=parameterfile)
#+++++#------  #
#+++++#------  # Check for Fortran90 Compiler
#+++++#------  #
#+++++#------  AC_CHECK_PROGS(COMPILER, $AbsoftF90 $CompaqF90 $IBMF90,no)
#+++++#------  if test "$COMPILER" = "no"; then 
#+++++#------     AC_MSG_ERROR(no valid compiler (f90,fort,xlf90) found)
#+++++#------  fi	
#+++++#------  
#+++++#------  #
#+++++#------  # define COMPILE and LINK VARIABLES
#+++++#------  #
#+++++#------  # Absoft f90
#+++++#------  if test "$COMPILER" = "$AbsoftF90"; then
#+++++#------     LINK="$COMPILER $absLFLAGS $USERLFLAGS"
#+++++#------     LLIBS="-lU77 -lm $USERLIBS"
#+++++#------     COMPILE="$COMPILER $absFLAGS $USERFLAGS"
#+++++#------     COMPILE77="$AbsoftF77 $absf77FLAGS"
#+++++#------     CPPVAR="$CPPVAR -DCPPVAR_U77 -DCPPVAR_ENDIANCHECK"
#+++++#------     ARCH=$absARCH
#+++++#------     FEXT=$absFEXT
#+++++#------     OPTIMIZE="$absOPTIMIZE90" 
#+++++#------     F77OPTIMIZE="$absOPTIMIZE77" 
#+++++#------     DEBUG="$absDEBUG90"      
#+++++#------     F77DEBUG="$absDEBUG77"      
#+++++#------  fi
#+++++#------  # Compaq fort
#+++++#------  if test "$COMPILER" = "$CompaqF90"; then
#+++++#------     LINK="$COMPILER $fortLFLAGS $USERLFLAGS"
#+++++#------     COMPILE="$COMPILER $fortFLAGS $USERFLAGS"
#+++++#------     LLIBS="$USERLIBS"
#+++++#------     COMPILE77="$CompaqF77 $fort77FLAGS"
#+++++#------     CPPVAR="$CPPVAR"
#+++++#------     ARCH=$fortARCH
#+++++#------     FEXT=$fortFEXT
#+++++#------     OPTIMIZE="$fortOPTIMIZE90"
#+++++#------     F77OPTIMIZE="$fortOPTIMIZE77"
#+++++#------     DEBUG="$fortDEBUG90"
#+++++#------     F77DEBUG="$fortDEBUG77"
#+++++#------  fi
#+++++#------  # IBM xlf90
#+++++#------  if test "$COMPILER" = "$IBMF90"; then
#+++++#------     LINK="$COMPILER $xlf90LFLAGS $USERLFLAGS"
#+++++#------     COMPILE="$COMPILER $xlf90FLAGS $USERFLAGS"
#+++++#------     COMPILE77="$IBMF77 $xlfFLAGS $USERFLAGS" 
#+++++#------     LLIBS="$USERLIBS"
#+++++#------     CPPVAR="$CPPVAR -DCPPVAR_SUPPORT_XLF -DCPPVAR_USAGE_EXIST -DCPPVAR_LANGEXT_XLF"
#+++++#------     ARCH=$xlf90ARCH
#+++++#------     FEXT=$xlf90FEXT
#+++++#------     OPTIMIZE="$xlfOPTIMIZE90"
#+++++#------     F77OPTIMIZE="$xlfOPTIMIZE77"
#+++++#------     DEBUG="$xlfDEBUG90"
#+++++#------     F77DEBUG="$xlfDEBUG77"
#+++++#------  fi
#+++++  #
#+++++  # setting architecture dependens directories
#+++++  #
#+++++  OBJDIR="$OBJDIR/$ARCH/$OBJECTS"
#+++++  PAWX="$BINDIR/$ARCH"
#+++++  MPICHROOT="$MPICHROOT-$ARCH"
#+++++  FFTWROOT="$FFTWROOT-$ARCH"
#+++++  #
#+++++  # use external libraries?
#+++++  #
#+++++  if test "$USELIBS" = "no"; then
#+++++    #
#+++++    # no external libraries according to command line argument
#+++++    #
#+++++    LLINK="$LLINK -L$OBJDIR"
#+++++    LLIBS="$LLIBS -lblas -llapack"
#+++++    CPPVAR="$CPPVAR -DCPPVAR_FFT_PACK"
#+++++    AC_MSG_RESULT(----------------------------------)
#+++++    AC_MSG_RESULT(according to command line argument)
#+++++    AC_MSG_RESULT(no external libraries are used    )
#+++++    AC_MSG_RESULT(----------------------------------)
#+++++      
#+++++#------  else
#+++++########################## beginning of absoft block###########################
#+++++#------    #
#+++++#------    # Determine external libraries for Absoft f90
#+++++#------    #
#+++++#------    if test "$COMPILER" = "$AbsoftF90"; then
#+++++#------      #
#+++++#------      # The FFTW-routines
#+++++#------      #
#+++++#------      # check for command line argument
#+++++#------      if test "$FFTWDIR" = "notset"; then
#+++++#------        # if FFTWDIR not set by command line argument
#+++++#------        # --> search default path (set with the variables)
#+++++#------        AC_MSG_RESULT(---------------------------------)
#+++++#------        AC_MSG_RESULT(take default path for FFTW as no )
#+++++#------        AC_MSG_RESULT(command line argument is given)
#+++++#------        AC_MSG_RESULT(---------------------------------)
#+++++#------        AC_CHECK_FILES($PAWDIR/$LIBRARIES/$FFTWROOT \
#+++++#------                       $PAWDIR/$LIBRARIES/$FFTWROOT/fftw/.libs/libfftw.a, 
#+++++#------                       FFTWDIR="$PAWDIR/$LIBRARIES/$FFTWROOT", 
#+++++#------                       FFTWDIR=notset)
#+++++#------      else
#+++++#------        # if FFTWDIR set by command line argument				   
#+++++#------        AC_CHECK_FILES($FFTWDIR/ \
#+++++#------                       $FFTWDIR/fftw/.libs/libfftw.a,
#+++++#------                       FFTWDIR="$FFTWDIR", 
#+++++#------       	               FFTWDIR=notset)
#+++++#------      fi			
#+++++#------      # if no FFTW-routines found - switch back to FFTPack
#+++++#------      # otherwise set variables
#+++++#------      if test "$FFTWDIR" = "notset"; then
#+++++#------        AC_MSG_WARN(-----------------------)
#+++++#------        AC_MSG_WARN(FFTW-library not found!)
#+++++#------        AC_MSG_WARN(using FFTPack instead  )
#+++++#------        AC_MSG_WARN(-----------------------)
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_FFT_PACK"
#+++++#------      else
#+++++#------        LLINK="$LLINK -L$FFTWDIR/fftw/.libs/"				     
#+++++#------        LLIBS="$LLIBS -lfftw"
#+++++#------        FFTW="yes"
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_FFT_FFTW"
#+++++#------      fi  
#+++++#------      #
#+++++#------      # The ATLAS BLAS routines
#+++++#------      #
#+++++#------      if test "$ATLASDIR" = "notset"; then
#+++++#------        # if ATLASDIR not set by command line argument
#+++++#------        if test "$ATLASARCH" = "notset"; then
#+++++#------          # if ATLASARCH not set by command line argument
#+++++#------          # the first in possible choice is taken 
#+++++#------          ATLASARCH=`ls $PAWDIR/$LIBRARIES/$ATLASROOT/lib/ \
#+++++#------                     | cut -d' ' -f1 | head -1`
#+++++#------          AC_MSG_RESULT(-------------------------------------------------)	   
#+++++#------          AC_MSG_RESULT(as no command line argument is given)
#+++++#------          AC_MSG_RESULT({$ATLASARCH} was chosen as architecture for ATLAS)	   
#+++++#------          AC_MSG_RESULT(if this is not the right choice - rerun )	   
#+++++#------          AC_MSG_RESULT(configure and use --with-atlasarch=... argument)	   
#+++++#------          AC_MSG_RESULT(-------------------------------------------------)	   
#+++++#------        fi
#+++++#------        AC_MSG_RESULT(---------------------------------)
#+++++#------        AC_MSG_RESULT(take default path for ATLAS as no)
#+++++#------        AC_MSG_RESULT(command line argument is given)
#+++++#------        AC_MSG_RESULT(---------------------------------)
#+++++#------        AC_CHECK_FILES($PAWDIR/$LIBRARIES/$ATLASROOT \
#+++++#------                       $PAWDIR/$LIBRARIES/$ATLASROOT/lib/$ATLASARCH/libf77blas.a\
#+++++#------                       $PAWDIR/$LIBRARIES/$ATLASROOT/lib/$ATLASARCH/libatlas.a, 
#+++++#------                       ATLASDIR="$PAWDIR/$LIBRARIES/$ATLASROOT", 
#+++++#------                       ATLASDIR=notset)
#+++++#------      else
#+++++#------        # if ATLASDIR set by command line argument				   
#+++++#------        if test "$ATLASARCH" = "notset"; then
#+++++#------          # if ATLASARCH not set by command line argument
#+++++#------          # the first in 
#+++++#------          ATLASARCH=`ls $PAWDIR/$LIBRARIES/$ATLASROOT/lib/ \
#+++++#------                    | cut -d' ' -f1 | head -1`
#+++++#------          AC_MSG_RESULT(-------------------------------------------------)	   
#+++++#------          AC_MSG_RESULT(as no command line argument is given)
#+++++#------          AC_MSG_RESULT({$ATLASARCH} was chosen as architecture for ATLAS)	   
#+++++#------          AC_MSG_RESULT(if this is not the right choice, rerun )	   
#+++++#------          AC_MSG_RESULT(reconfigure and use --with-atlasarch=... argument)	   
#+++++#------          AC_MSG_RESULT(-------------------------------------------------)	   
#+++++#------        fi
#+++++#------        AC_CHECK_FILES($ATLASDIR/ \
#+++++#------                      $ATLASDIR/lib/$ATLASARCH/libf77blas.a\
#+++++#------                      $ATLASDIR/lib/$ATLASARCH/libatlas.a,
#+++++#------                      ATLASDIR="$ATLASDIR", 
#+++++#------                      ATLASDIR=notset)
#+++++#------      fi			
#+++++#------      # if no ATLAS-library found - switch back to internal blas
#+++++#------      # otherwise set variables
#+++++#------      if test "$ATLASDIR" = "notset"; then
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        AC_MSG_WARN(ATLAS-library not found!)
#+++++#------        AC_MSG_WARN(using internal blas instead)
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        LLIBS="$LLIBS -llapack -lblas"
#+++++#------      else
#+++++#------        LLINK="$LLINK -L$ATLASDIR/lib/$ATLASARCH"				     
#+++++#------        LLIBS="$LLIBS -llapack -lf77blas -latlas"
#+++++#------      fi  
#+++++#------    fi  
#+++++########################## end of absoft block###########################
#+++++
#+++++########################## beginning of ibm block###########################
#+++++#------    # 
#+++++#------    # Determine external libraries for IBM xlf90
#+++++#------    #
#+++++#------    if test  "$COMPILER" = "$IBMF90"; then
#+++++#------      # THE ESSL 
#+++++#------      AC_HAVE_LIBRARY(essl,ESSL=found,ESSL=no)
#+++++#------      if test "$ESSL" = "no"; then
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        AC_MSG_WARN(ESSL library not found)
#+++++#------        AC_MSG_WARN(using internal blas instead)
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        LLIBS="$LLIBS -lblas -llapack"
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_FFT_PACK"
#+++++#------        
#+++++#------                              else
#+++++#------        LLIBS="$LLIBS -lessl"			    
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_BLAS_ESSL -DCPPVAR_FFT_ESSL"
#+++++#------      fi   
#+++++#------    fi
#+++++########################## end of ibm block###########################
#+++++  
#+++++########################## beginning of Compaq###########################
#+++++#------    # 
#+++++#------    # Determine external library for Compaq fort
#+++++#------    #
#+++++#------    if test  "$COMPILER" = "$CompaqF90"; then
#+++++#------      # THE CXML
#+++++#------      AC_HAVE_LIBRARY(cxml,CXML=found,CXML=no)
#+++++#------      if test "$CXML" = "no"; then
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        AC_MSG_WARN(CXML library not found)
#+++++#------        AC_MSG_WARN(using internal blas instead)
#+++++#------        AC_MSG_WARN(---------------------------)
#+++++#------        LLIBS="$LLIBS -lblas -llapack"
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_FFT_PACK"
#+++++#------                              else
#+++++#------        LLIBS="$LLIBS -lcxml"	 		    
#+++++#------        CPPVAR="$CPPVAR -DCPPVAR_FFT_PACK"
#+++++#------      fi   
#+++++#------    fi
#+++++########################## end of ibm block###########################
#+++++  fi
#+++++
#+++++  LLIBS="$LLIBS $USERLIBS"
#+++++  #
#+++++  # Check for MPICH
#+++++  #
#+++++  # only if non-IBM!!!
#+++++  #
#+++++#------  if test "$ARCH" = "$xlf90ARCH"; then
#+++++#------    # assume, that libmpi.a is in search path or has 
#+++++#------    #AC_HAVE_LIBRARY(mpi,LLIBSPARALLEL="$LLIBS -lmpi",LLIBSPARALLEL="$LLIBS")
#+++++#------    LLIBSPARALLEL="$LLIBS -lmpi"
#+++++#------  else
#+++++    if test "$MPICHDIR" = "notset"; then
#+++++      # if MPICHDIR not set by command line argument
#+++++      AC_MSG_RESULT(---------------------------------)
#+++++      AC_MSG_RESULT(take default path for MPICH as no)
#+++++      AC_MSG_RESULT(command line argument is given)
#+++++      AC_MSG_RESULT(---------------------------------)
#+++++      AC_CHECK_FILES($PAWDIR/$LIBRARIES/$MPICHROOT \
#+++++                     $PAWDIR/$LIBRARIES/$MPICHROOT/lib/libmpich.a\
#+++++                     $PAWDIR/$LIBRARIES/$MPICHROOT/lib/libfmpich.a,
#+++++                     MPICHDIR="$PAWDIR/$LIBRARIES/$MPICHROOT", 
#+++++                     MPICHDIR=notset)
#+++++    else
#+++++      # if MPICHDIR set by command line argument				   
#+++++      AC_CHECK_FILES($MPICHDIR/ \
#+++++                     $MPICHDIR/lib/libmpich.a\
#+++++                     $MPICHDIR/lib/libfmpich.a,
#+++++                     MPICHDIR="$MPICHDIR", 
#+++++                     MPICHDIR=notset)
#+++++    fi			
#+++++    # if no MPICH-library found -> no parallel compilation possible
#+++++    # otherwise set variables
#+++++    if test "$MPICHDIR" = "notset"; then
#+++++      AC_MSG_WARN(--------------------------------)
#+++++      AC_MSG_WARN(MPICH-library not found!)
#+++++      AC_MSG_WARN(no parallel compilation possible)
#+++++      AC_MSG_WARN(--------------------------------)
#+++++    else
#+++++      LLINKPARALLEL="$LLINK -L$MPICHDIR/lib"				     
#+++++      LLIBSPARALLEL="$LLIBS -lmpich -lfmpich"
#+++++    fi  
#------  fi
fi     # end of block "without parmfile" #--------------------------------------

LLIBSPARALLEL="$LLIBSPARALLEL $USERLIBS"
#
# introduction of PARALLELNAME to avoid problems
# with Makefile.* generation when PARALLELDIR=fast
# --> see below 
if test "$PARALLELDIR" = "$FASTDIR"; then
  PARALLELNAME="parallel"
else
  PARALLELNAME="$PARALLELDIR"
fi
#
# Definition of variables to be substituted in the *.in files
#   There the synax is @VARIABLE@
#
#ATTENTION CHANGE!!!!
CPP="$CPP -traditional"

#
#########################################################################
##                                                                     ##
##  define substitutions for AC_OUTPUT                                 ##
##                                                                     ##
#########################################################################
AC_SUBST(PAWDIR)
AC_SUBST(COMPILE)
AC_SUBST(COMPILE77)
AC_SUBST(LINK)
AC_SUBST(LLINK)
AC_SUBST(LLINKPARALLEL)
AC_SUBST(LLIBSPARALLEL)
AC_SUBST(LLIBS)
AC_SUBST(OBJDIR)
AC_SUBST(SRCDIR)
AC_SUBST(BINDIR)
AC_SUBST(FFTW)
AC_SUBST(FFTWDIR)
AC_SUBST(LAPACKDIR)
AC_SUBST(MPICHDIR)
AC_SUBST(TOOLDIR)
AC_SUBST(NONEDIR)
AC_SUBST(FASTDIR)
AC_SUBST(DEBUGDIR)
AC_SUBST(PROFDIR)
AC_SUBST(PARALLELDIR)
AC_SUBST(PARALLELNAME)
AC_SUBST(TOOLCOMP)
AC_SUBST(PAWX)
AC_SUBST(CPPVAR)
AC_SUBST(FEXT)
AC_SUBST(ARCH)
AC_SUBST(CPP)
AC_SUBST(MAKE)
AC_SUBST(xlf90ARCH)
#
#########################################################################
##                                                                     ##
##  create Makefiles and f90pp                                         ##
##                                                                     ##
#########################################################################
#
# create Makefile from Makefile.in 
#
#C_OUTPUT(Makefile:Makefile.in $BINDIR/f90pp:$SOURCES/F90PP/f90pp.in)
AC_OUTPUT(Makefile:Makefile.in)
#
# create f90pp from f90pp.in and make it executable
#
AC_OUTPUT($BINDIR/f90pp:$PAWDIR/src/F90PP/f90pp.in)
chmod +x $BINDIR/f90pp
#
# substituting to get Makefile.none
#
AC_MSG_RESULT(creating Makefile.none)
sed -e "s%\(OBJDIR.*=.*\)%\1/$NONEDIR%" \
    -e "s%\(PAWX.*=.*\)%\1/$NONEBIN%"\
    -e "s%varOPTRM%$NONEOPTRM%"\
    -e 's%varPARALLEL%F%'\
    Makefile > Makefile.${NONEDIR}
#
# subsituting to get Makefile.fast
#
AC_MSG_RESULT(creating Makefile.fast)
sed -e "s%\(COMPILE=.*\)%\1 $OPTIMIZE %"\
    -e "s%\(COMPILE77.*=.*\)%\1 $F77OPTIMIZE %"\
    -e "s%\(OBJDIR.*=.*\)%\1/$FASTDIR%"\
    -e "s%\(PAWX.*=.*\)%\1/$FASTBIN%"\
    -e "s%varOPTRM%$FASTOPTRM%"\
    -e 's%varPARALLEL%F%'\
    Makefile > Makefile.${FASTDIR}
#
# subsituting to get Makefile.parallel
#
AC_MSG_RESULT(creating Makefile.parallel)
sed -e "s%\(COMPILE=.*\)%\1 $OPTIMIZE %"\
      -e "s%\(COMPILE77.*=.*\)%\1 $F77OPTIMIZE %"\
      -e "s%\(OBJDIR.*=.*\)%\1/$PARALLELDIR%"\
      -e "s%\(PAWX.*=.*\)%\1/$PARALLELBIN%"\
      -e "s%varOPTRM%$PARALLELOPTRM%"\
      -e 's%varPARALLEL%T%'\
      Makefile > Makefile.${PARALLELNAME}
#
# subsituting to get Makefile.debug
#
AC_MSG_RESULT(creating Makefile.debug)
sed -e "s%\(COMPILE=.*\)%\1 $DEBUG %"\
    -e "s%\(COMPILE77=.*\)%\1 $F77DEBUG %"\
    -e "s%\(OBJDIR.*=.*\)%\1/$DEBUGDIR%"\
    -e "s%\(PAWX.*=.*\)%\1/$DEBUGBIN%"\
    -e "s%varOPTRM%$DEBUGOPTRM%"\
    -e 's%varPARALLEL%F%'\
    Makefile > Makefile.${DEBUGDIR}
#
# subsituting to get Makefile.prof
#
AC_MSG_RESULT(creating Makefile.prof)
sed -e "s%\(COMPILE=.*\)%\1 $PROF %"\
    -e "s%\(COMPILE77=.*\)%\1 $F77PROF %"\
    -e "s%\(LINK=.*\)%\1 -pg %"\
    -e "s%\(OBJDIR.*=.*\)%\1/$PROFDIR%"\
    -e "s%\(PAWX.*=.*\)%\1/$PROFBIN%"\
    -e "s%varOPTRM%$PROFOPTRM%"\
    -e 's%varPARALLEL%F%'\
    Makefile > Makefile.${PROFDIR}
#
# creating subdirectories
#
AC_MSG_RESULT(creating subdirectories and copying shell scripts)
mkdir -p  ${BINDIR} ${BINDIR}/${ARCH} 
list="$NONEDIR $FASTDIR $DEBUGDIR $PARALLELDIR $PROFDIR"
for i in $list; do
  mkdir -p ${OBJDIR}/${i}
done
#
#  copy script files from
#
cp $SCRIPTDIR/* $BINDIR/$ARCH/
#
# copy current parameter file to "parms.in_use"
#
AC_MSG_RESULT(copy $PARMFILE to parms.in_use)
cp $PARMFILE parms.in_use
#
# removal of config.cache 
#
rm -f config.cache
#########################################################################
##                                                                     ##
##                      print summary                                  ##
##                                                                     ##
#########################################################################
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(printing out of summary              )
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(architecture                   : $ARCH)
AC_MSG_RESULT(CPP command                    : $CPP)
AC_MSG_RESULT(MAKE command                   : $MAKE)
AC_MSG_RESULT(F90 compile command            : $COMPILE)
AC_MSG_RESULT(F77 compile command            : $COMPILE77)
AC_MSG_RESULT(F90 optimize options           : $OPTIMIZE)
AC_MSG_RESULT(F77 optimize options           : $F77OPTIMIZE)
AC_MSG_RESULT(F90 debug flags                : $DEBUG)
AC_MSG_RESULT(F77 debug flags                : $F77DEBUG)
AC_MSG_RESULT(F90 profiling flags            : $PROF)
AC_MSG_RESULT(F77 profiling flags            : $F77PROF)
AC_MSG_RESULT(F90 linking flags              : $LFLAGS)
AC_MSG_RESULT(search directories (sequential): $LLINK)
AC_MSG_RESULT(search directories (parallel)  : $LLINKPARALLEL)
AC_MSG_RESULT(external libraries (sequential): $LLIBS)
AC_MSG_RESULT(external libraries (parallel)  : $LLIBSPARALLEL)
AC_MSG_RESULT(preprocessor variables         : $CPPVAR)
AC_MSG_RESULT(F90 file extension             : $FEXT)
AC_MSG_RESULT(lapack/blas library            : $LAPACKDIR)
AC_MSG_RESULT(libs for Fourier transforms    : $FFTWDIR)
AC_MSG_RESULT(mpi library                    : $FFTWDIR)
AC_MSG_RESULT(-------------------------------------)
AC_MSG_RESULT(configuration completed successfully!)


































































