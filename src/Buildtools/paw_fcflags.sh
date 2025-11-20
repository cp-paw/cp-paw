#/bin/bash
export THISDIR=$(pwd)
#________1_________2_________3_________4_________5_________6_________7_________8
################################################################################
##     
##    paw_fcflags.sh   
##     
##   returns a set of compile flags from
##
## https://fortran-lang.discourse.group/t/compilation-flags-advice-for-production-and-distribution/2821/2
##     
##   Remark: compile option -c is not included to permit compiling and linking
##           in one step.  
##     
##   Author P. Bloechl May, 2024 
################################################################################
# 
#-------------------------------------------------------------------------------
#  help message
#-------------------------------------------------------------------------------
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \t paw_fcflags.sh options \n"
USAGE="$USAGE \n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -s (Darwin,Linux, Windows) operating system\n"
USAGE="$USAGE \t -c (gfortran, ifort) Fortran compiler\n"
USAGE="$USAGE \t -t (release, debug) type of compilation flags \n"
USAGE="$USAGE \t -v verbose (false) causes error exit \n"
USAGE="$USAGE \t -h prints this help message\n"
USAGE="$USAGE \n"
USAGE="$USAGE -c is the suffix if no parmfile is present or -f ''\n"

#-------------------------------------------------------------------------------
#  Resolve arguments
#-------------------------------------------------------------------------------
export OS=
export COMPILER=
export TYPE=
export VERBOSE=false
while getopts :s:c:t:vh OPT ; do
  case $OPT in
    s) OS=$OPTARG ;;
    c) COMPILER=$OPTARG ;;
    t) TYPE=$OPTARG ;;
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
if [[ -z $OS ]] ; then OS=$(uname -s) ; fi
if [[ -z $COMPILER ]] ; then 
  echo "error in $0: no compiler specified option -c is mandatory"           >&2
  exit 1
else
  if [[ -z $(echo "ifort gfortran" | grep ${COMPILER}) ]]; then
    echo "error in $0: invalid compiler selection"                           >&2
    echo "allowed are ifort gfortran"                                        >&2
    exit 1
  fi
fi
if [[ -z $TYPE ]] ; then 
  echo "error in $0: no type specified. option -t is mandatory"              >&2
  exit 1
else
  if [[ -z $(echo "debug release" | grep ${TYPE}) ]]; then
    echo "error in $0: invalid type selected with -t"                        >&2
    echo "allowed are debug release"                                         >&2
    exit 1
  fi
fi
if [[ $VERBOSE = true ]] ; then
  echo "Caution: option -v is only for testing and will cause an error exit"
  echo "Operating system = $OS" 
  echo "Compiler         = $COMPILER" 
  echo "Type             = $TYPE" 
fi

# I don’t remember the specific reasons, but I ended up separating the
# GNU release flags for macOS from Linux/Windows because a particular
# unknown flag, switched on by -O3 (or -flto), led to segfaults on
# MacOS. This happened several years ago with GNU 7/8/9(?). The bugs
# may have been resolved in the newer releases of GNU compilers.

# These flags are taken from the CMake files of the ParaMonte library
# 5. These flags do not include architecture-specific optimization
# flags that Steve Kargl listed (to ensure portability of the
# generated library).

#Perhaps, a compilation of the flags similar to the above should
#appear on the FortranLang website, if there is not any there already.

# Intel has an excellent summary 9 of its compiler flags.

# One final note, enabling interprocedural optimizations (e.g., with
# -flto or -ipo) will significantly lengthen the compilation process
# possibly extending a 30-sec process to 30 mins.

# I mentioned only gfortran and ifort because I have experience
# primarily with these two.


#  OS can be Unix, Darwin, Windows
#  COMPILER can be ifort, gfortran
#  TYPE can be debug, release
#
#  Unix          ifort debug
#  Windows       ifort debug
#  Windows       ifort release
#  Unix          ifort release
#   *            gfortran debug
#  Darwin        gfortran release
#  Linux/Windows gfortran release
#
#

if [[ $OS = Unix && $COMPILER = ifort && $TYPE=debug ]] ; then
  #-debug full           # generate full debug information
  #-g3                   # generate full debug information
  #-O0                   # disable optimizations
  #-CB                   # Perform run-time bound-checks on array subscript 
  #                      # and substring references (same as the -check bounds 
  #                      $ option)
  #-init:snan,arrays     # initialize arrays and scalars to NaN
  #-warn all             # enable all warning
  #-gen-interfaces       # generate interface block for each routine 
  #                      # in the source files
  #-traceback            # trace back for debugging
  #-check all            # check all
  #-check bounds         # check array bounds
  #-fpe-all=0            # Floating-point invalid, divide-by-zero, 
  #                      # and overflow exceptions are enabled
  #-fpe0                 # Ignore underflow (yield 0.0); Abort on other 
  #                        IEEE exceptions.
  #-diag-error-limit=10  # max diagnostic error count
  #-diag-disable=5268    # Extension to standard: The text exceeds right hand 
  #                        column allowed on the line.
  #-diag-disable=7025    # This directive is not standard Fxx.
  #-diag-disable=10346   # optimization reporting will be enabled at link time 
  #                        when performing interprocedural optimizations.
  #-ftrapuv             # Initializes stack local variables to an unusual value 
  #                        to aid error detection.
  # FCFLAGS="-debug full -g3 -O0 -CB -init:snan,arrays -warn all \
  #          -gen-interfaces -traceback \
  #          -check all -check bounds \
  #          -fpe-all=0 -fpe0 \
  #          -diag-error-limit=10 -diag-disable=5268 -diag-disable=7025 \
  #          -diag-disable=10346 -ftrapuv "

  # flags used by axel ehrich 
  # FCFLAGS="-c -O0 -CA -CU -CB -g \
  #          -check udio_iostat -check stack 
  #          -check output_conversion -check contiguous -check assume \
  #          -check bounds -check arg_temp_created -check format \
  #          -check uninit -check bounds -check format \
  #          -check pointers -check uninit \
  #          -debug full -debug-parameters all \
  #          -fp-stack-check -traceback -warn declarations"
  # compromise between axel and default
  FCFLAGS="-debug full -g2 -O0 -CB -init:snan,arrays -warn all \
           -gen-interfaces -traceback \
           -check udio_iostat -check stack 
           -check output_conversion -check contiguous -check assume \
           -check bounds -check arg_temp_created -check format \
           -check uninit -check bounds -check format \
           -fpe-all=0 -fpe0 \
           -diag-error-limit=10 -diag-disable=5268 -diag-disable=7025 \
           -diag-disable=10346 -ftrapuv "

elif [[ $OS = Windows && $COMPILER = ifort && $TYPE = debug ]] ; then
  # /debug:full
  # /stand:f18      # issue compile-time messages for nonstandard 
  #                   language elements.
  # /Zi
  # /CB
  # /Od
  # /Qinit:snan,arrays
  # /warn:all
  # /gen-interfaces
  # /traceback
  # /check:all
  # /check:bounds
  # #/fpe-all:0
  # /fpe:0
  # /Qdiag-error-limit:10
  # /Qdiag-disable:5268
  # /Qdiag-disable:7025
  # /Qtrapuv
  #
  FCFLAGS="/debug:full /stand:f18 /Zi /CB /Od /Qinit:snan,arrays /warn:all \
           /gen-interfaces /traceback /check:all /check:bounds  /fpe:0 \
           /Qdiag-error-limit:10 /Qdiag-disable:5268 /Qdiag-disable:7025 \
           /Qtrapuv"

           #/fpe-all:0

 elif [[ $OS = Windows && $COMPILER = ifort && $TYPE = release ]] ; then
  # /O3                     # Enable O3 optimization.
  # /Qvec                   # enable vectorization.
  # /Qunroll                # [:n] set the maximum number of times 
  #                           to unroll loops (no number n means automatic).
  # /Qunroll-aggressive     # use more aggressive unrolling for certain loops.
  # /Qinline-forceinline    # Instructs the compiler to force inlining of 
  #                           functions suggested for inlining whenever 
  #                           the compiler is capable doing so.
  # #/Qguide-vec:4          # enable guidance for auto-vectorization, 
  #                           causing the compiler to generate messages 
  #                           suggesting ways to improve optimization 
  #                           (default=4, highest).
  # #/Qparallel             # generate multithreaded code for loops that 
  #                           can be safely executed in parallel.
  # #/Qipo-c:               # Tells the compiler to optimize across multiple 
  #                           files and generate a single object file 
  #                           ipo_out.obj without linking
  #                           info at: https://software.intel.com/en-us/Fortran-compiler-developer-guide-and-reference-ipo-c-qipo-c
  # /Qftz                   #
  # /Qipo                   # enable interprocedural optimization between files.
  # /Qip                    # 

  FCFLAGS="/O3 /Qvec /Qunroll /Qunroll-aggressive /Qinline-forceinline 
          /Qftz /Qipo /Qip"

elif [[ $OS = Unix && $COMPILER = ifort && $TYPE = release ]] ; then
  #  -stand f18                # issue compile-time messages for 
  #                              nonstandard language elements.
  #  -O3                       # set the optimizations level
  #  -unroll                   # [=n] set the maximum number of times to 
  #                              unroll loops (no number n means automatic).
  #  -unroll-aggressive        # use more aggressive unrolling for certain 
  #                              loops.
  #  -diag-disable=10346       # optimization reporting will be enabled at 
  #                              link time when performing interprocedural 
  #                              optimizations.
  #  -diag-disable=10397       # optimization reporting will be enabled at 
  #                              link time when performing interprocedural 
  #                              optimizations.
  #  #-guide-vec=4             # enable guidance for auto-vectorization, 
  #                              causing the compiler to generate messages 
  #                              suggesting ways to improve optimization 
  #                              (default=4, highest).
  #  #-parallel                # generate multithreaded code for loops that 
  #                              can be safely executed in parallel. This 
  #                              option requires MKL libraries.
  #  #-qopt-subscript-in-range # assumes there are no "large" integers being 
  #                              used or being computed inside loops. 
  #                              A "large" integer is typically > 2^31.
  #  -ftz                      # Flushes subnormal results to zero.
  #  -inline-forceinline       # Instructs the compiler to force inlining of 
  #                              functions suggested for inlining whenever 
  #                              the compiler is capable doing so.
  #  -finline-functions        # enables function inlining for single file 
  #                              compilation.
  #  -ipo                      # enable interprocedural optimization 
  #                            between files.
  #   -ip        
  
  FCFLAGS="-stand f18 -O3 -unroll -unroll-aggressive \
           -diag-disable=10346  -diag-disable=10397 \
           -ftz -inline-forceinline -finline-functions-ipo -ip"

  # axels parameter
  # FCFLAGS="-c -O3 -no-ipo -no-ip -xCORE-AVX2 -finline-functions \
  #           -finline-limit=50"

elif [[ $COMPILER = gfortran && $TYPE = debug ]] ; then
  #  -g3                         # generate full debug information
  #  -O0                         # disable optimizations
  # #-fsanitize=undefined        # enable UndefinedBehaviorSanitizer 
  #                                for undefined behavior detection.
  # #-fsanitize=address          # enable AddressSanitizer, for memory 
  #                                error detection, like out-of-bounds 
  #                                and use-after-free bugs.
  # #-fsanitize=leak             # enable LeakSanitizer for memory leak detection.
  #  -fcheck=all                 # enable the generation of run-time checks
  #  -ffpe-trap=invalid,zero,overflow  # ,underflow : Floating-point invalid,
  #                                      divide-by-zero, and overflow 
  #                                      exceptions are enabled
  #  -ffpe-summary=all           # Specify a list of floating-point exceptions,
  #                                whose flag status is printed to ERROR_UNIT
  #                                when invoking STOP and ERROR STOP.
  #                              # Can be either 'none', 'all' or a  $
  #                                comma-separated list of the following 
  #                                exceptions: 'invalid','zero','overflow',
  #                                'underflow', 'inexact' and 'denormal'
  #  -finit-integer=-2147483647  # initilize all integers to negative infinity
  #  -finit-real=snan            # initialize REAL and COMPLEX variables 
  #                                with a signaling NaN
  #  -fbacktrace                 # trace back for debugging
  # #-pedantic                   # issue warnings for uses of extensions to 
  #                                the Fortran standard. Gfortran10 with 
  #                                MPICH 3.2 in debug mode crashes with this 
  #                                flag at mpi_bcast. Excluded until MPICH #
  #                                upgraded.
  #  -fmax-errors=10             # max diagnostic error count
  #  -Wno-maybe-uninitialized    # avoid warning of no array pre-allocation.
  #  -Wall                       # enable all warnings:
  #                              # -Waliasing, -Wampersand, -Wconversion, 
  #                                -Wsurprising, -Wc-binding-type, 
  #                                -Wintrinsics-std, -Wtabs, 
  #                                -Wintrinsic-shadow, -Wline-truncation, 
  #                                -Wtarget-lifetime, -Winteger-division, 
  #                                -Wreal-q-constant, -Wunused, 
  #                                -Wundefined-do-loop
  #                              # gfortran10 crashes and cannot compile MPI 
  #                                ParaMonte with mpich in debug mode. 
  #                                Therefore -wall is disabled for now, until
  #                                 MPICH upgrades interface.
  # #-Wconversion-extra          # Warn about implicit conversions between 
  #                                different types and kinds. This option 
  #                                does not imply -Wconversion.
  # #-Werror=conversion          # Turn all implicit conversions into an 
  #                                error. This is important to avoid 
  #                                inadvertent implicit change of precision 
  #                                in generic procedures of various kinds, 
  #                                due to the use of `RK` to represent 
  #                                different kinds.
  # #-Werror=conversion-extra    # Turn all implicit conversions into an error.
  #                                This is too aggressive and as such currently
  #                                deactivated. For example, it yields an error
  #                                on the multiplication of integer with real.
  #  -fno-unsafe-math-optimizations #
  #  -fsignaling-nans            #
  #  -frounding-math             #
  #  -Wno-surprising             #

  FCFLAGS="-g3 -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow \
           -ffpe-summary=all \
           -finit-integer=-2147483647  -finit-real=snan -fbacktrace \
           -fmax-errors=10 -Wno-maybe-uninitialized -Wall \
           -fno-unsafe-math-optimizations \
           -fsignaling-nans -frounding-math -Wno-surprising \
           -march=native"


elif [[ $OS = Darwin && $COMPILER = gfortran && $TYPE = release ]] ; then
  # -fauto-inc-dec
  # -fbranch-count-reg
  # -fcombine-stack-adjustments
  # -fcompare-elim
  # -fcprop-registers
  # -fdce
  # -fdefer-pop
  # - #-fdelayed-branch
  # -fdse
  # -fforward-propagate
  # -fguess-branch-probability
  # -fif-conversion
  # -fif-conversion2
  # -finline-functions-called-once
  # -fipa-profile
  # -fipa-pure-const
  # -fipa-reference
  # -fipa-reference-addressable
  # -fmerge-constants
  # -fmove-loop-invariants
  # -fomit-frame-pointer
  # -freorder-blocks
  # -fshrink-wrap
  # -fshrink-wrap-separate
  # -fsplit-wide-types
  # -fssa-backprop
  # -fssa-phiopt
  # -ftree-bit-ccp
  # -ftree-ccp
  # -ftree-ch
  # -ftree-coalesce-vars
  # -ftree-copy-prop
  # -ftree-dce
  # -ftree-dominator-opts
  # -ftree-dse
  # -ftree-forwprop
  # -ftree-fre
  # -ftree-phiprop
  # -ftree-pta
  # -ftree-scev-cprop
  # -ftree-sink
  # -ftree-slsr
  # -ftree-sra
  # -ftree-ter
  # -funit-at-a-time
  # -falign-functions   # -falign-jumps
  # -falign-labels   # -falign-loops
  # -fcaller-saves
  # -fcode-hoisting
  # -fcrossjumping
  # -fcse-follow-jumps   # -fcse-skip-blocks
  # -fdelete-null-pointer-checks
  # -fdevirtualize   # -fdevirtualize-speculatively
  # -fexpensive-optimizations
  # -fgcse   # -fgcse-lm
  # -fhoist-adjacent-loads
  # -finline-functions
  # -finline-small-functions
  # -findirect-inlining
  # -fipa-bit-cp   # -fipa-cp   # -fipa-icf
  # -fipa-ra   # -fipa-sra   # -fipa-vrp
  # -fisolate-erroneous-paths-dereference
  # -flra-remat
  # -foptimize-sibling-calls
  # -foptimize-strlen
  # -fpartial-inlining
  # -fpeephole2
  # -freorder-blocks-algorithm=stc
  # -freorder-blocks-and-partition   # -freorder-functions
  # -frerun-cse-after-loop
  # -fschedule-insns   # -fschedule-insns2
  # -fsched-interblock   # -fsched-spec
  # -fstore-merging
  # -fstrict-aliasing
  # -fthread-jumps
  # -ftree-builtin-call-dce
  # -ftree-pre
  # -ftree-switch-conversion   # -ftree-tail-merge
  # -ftree-vrp
  # -fgcse-after-reload
  # -fipa-cp-clone
  # -floop-interchange
  # -floop-unroll-and-jam
  # -fpeel-loops
  # -fpredictive-commoning
  # -fsplit-paths
  # -ftree-loop-distribute-patterns
  # -ftree-loop-distribution
  # -ftree-loop-vectorize
  # -ftree-partial-pre
  # -ftree-slp-vectorize
  # -funswitch-loops
  # -fvect-cost-model
  # -fversion-loops-for-strides

  FCFLAGS="-fauto-inc-dec -fbranch-count-reg \
 -fcombine-stack-adjustments -fcompare-elim -fcprop-registers -fdce \
 -fdefer-pop -fdse -fforward-propagate -fguess-branch-probability \
 -fif-conversion -fif-conversion2 -finline-functions-called-once \
 -fipa-profile -fipa-pure-const -fipa-reference \
 -fipa-reference-addressable -fmerge-constants -fmove-loop-invariants \
 -fomit-frame-pointer -freorder-blocks -fshrink-wrap \
 -fshrink-wrap-separate -fsplit-wide-types -fssa-backprop -fssa-phiopt \
 -ftree-bit-ccp -ftree-ccp -ftree-ch -ftree-coalesce-vars \
 -ftree-copy-prop -ftree-dce -ftree-dominator-opts -ftree-dse \
 -ftree-forwprop -ftree-fre -ftree-phiprop -ftree-pta \
 -ftree-scev-cprop -ftree-sink -ftree-slsr -ftree-sra -ftree-ter \
 -funit-at-a-time -falign-functions -falign-jumps -falign-labels \
 -falign-loops -fcaller-saves -fcode-hoisting -fcrossjumping \
 -fcse-follow-jumps -fcse-skip-blocks -fdelete-null-pointer-checks \
 -fdevirtualize -fdevirtualize-speculatively -fexpensive-optimizations \
 -fgcse -fgcse-lm -fhoist-adjacent-loads -finline-functions \
 -finline-small-functions -findirect-inlining -fipa-bit-cp -fipa-cp \
 -fipa-icf -fipa-ra -fipa-sra -fipa-vrp \
 -fisolate-erroneous-paths-dereference -flra-remat \
 -foptimize-sibling-calls -foptimize-strlen -fpartial-inlining \
 -fpeephole2 -freorder-blocks-algorithm=stc \
 -freorder-blocks-and-partition -freorder-functions \
 -frerun-cse-after-loop -fschedule-insns -fschedule-insns2 \
 -fsched-interblock -fsched-spec -fstore-merging -fstrict-aliasing \
 -fthread-jumps -ftree-builtin-call-dce -ftree-pre \
 -ftree-switch-conversion -ftree-tail-merge -ftree-vrp \
 -fgcse-after-reload -fipa-cp-clone -floop-interchange \
 -floop-unroll-and-jam -fpeel-loops -fpredictive-commoning \
 -fsplit-paths -ftree-loop-distribute-patterns \
 -ftree-loop-distribution -ftree-loop-vectorize -ftree-partial-pre \
 -ftree-slp-vectorize -funswitch-loops -fvect-cost-model \
 -fversion-loops-for-strides -march=native"

elif [[ ( $OS = Linux || $OS = Windows ) \
           && $COMPILER = gfortran && $TYPE = release ]] ; then
  # -ftree-vectorize        # perform vectorization on trees. enables  
  #                           -ftree-loop-vectorize and -ftree-slp-vectorize.
  # -funroll-loops          # [=n] set the maximum number of times to unroll 
  #                           loops (no number n means automatic).
  # -O3                     # set the optimizations level
  # -finline-functions      # consider all functions for inlining, even if 
  #                           they are not declared inline.
  # #-fwhole-program        # allow the compiler to make assumptions on the 
  #                           visibility of the symbols leading to more 
  #                           aggressive optimization decisions.
  # -flto=3                 #

 FCFLAGS="-ftree-vectorize -funroll-loops -O3 -finline-functions \
          -fwhole-program -flto=3 -march=native" 
fi

# because flags are written to sttout, VERBOSE messes up the result
if [[ $VERBOSE = true ]] ; then
  echo "intended error exit from $0: option -v only for debugging"           >&2
  echo FCFLAGS=$FCFLAGS                                                      >&2
  exit 1
fi
echo ${FCFLAGS} > /dev/stdout
