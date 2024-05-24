## #/bin/sh  ! comment this out for openmpi
################################################################################
################################################################################
####    parameters to be set by the user:
####
##     default for the number of nodes to be used. overwritten with option -n
export NNODES=10
##
##     paw executable for parallel computatiuon. Can be left as is if found
export PARCODE=$(which ppaw_fast.x)
##
##     path of the mpirun command. Can be left as is if found
export MPIRUN=$(which mpirun)          #package mpi execution
##
####  
################################################################################
################################################################################
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \tdoppaw(.sh) options JOB\n"
USAGE="$USAGE \n"
USAGE="$USAGE An ending ".cntl" of JOB is stripped away.\n"
USAGE="$USAGE \n"
USAGE="$USAGE Purpose:\n"
USAGE="$USAGE \t perform a parallel paw_calculation\n"
USAGE="$USAGE \t\t with control file JOB.cntl\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -n nnodes: submit parallel job with the specified number \
                             of processes\n"
USAGE="$USAGE \t -x name: name of executable ppaw_XXX.x with full path\n"
USAGE="$USAGE \t -h \t print this help message \n"
USAGE="$USAGE \t -0: dry-run (creates files but does not run jobs)\n"
USAGE="$USAGE \t -v: verbose\n"
USAGE="$USAGE \n"
USAGE="$USAGE attention: contains hardcoded user-specified parameters \n"
USAGE="$USAGE attention: creates and removes a temporary directory in \tmp \n"
USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
#  implement dry-run
#-------------------------------------------------------------------------------
function execute(){
    # || is "or" in [[...]], && is "and" in [[..]] 
    if [[ "${DRYRUN}" = "yes" || "${VERBOSE}" = "yes" ]] ; then
      echo "${@}"
    fi
    if [[ "${DRYRUN}" = "no" ]] ; then
      eval "$@"
    fi
}
#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
export DRYRUN=no
export VERBOSE=no
while getopts :hv0n:x: OPT  ; do
  case $OPT in
    n)   #nodes:
      NNODES=$OPTARG
      ;;
    x)   #name of the executable ppaw_fast.x with full path
      PARCODE=$OPTARG
      ;;
    0)   #nothing:
      DRYRUN=yes
      ;;
    v)   #verbose
      VERBOSE=yes
      ;;
    h)   # help
      echo -e $USAGE
      exit 0
      ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0" >&2
      echo "invalid option -$OPTARG" >&2
      echo "retrieve argument list with:" >&2
      echo "$0 -h" >&2
      exit 1
      ;;
    :)    # no argument passed to option requiring one
      echo "error in $0" >&2
      echo "option -$OPTARG requires an additional argument" >&2
      exit 1
      ;;  
  esac
done
shift $(($OPTIND - 1))
JOB=$1  #  Directory containing all subprojects relative to THISDIR
if [ -z "$JOB" ] ; then 
  echo "error in $0: project name missing"
  echo -e "$USAGE"
  exit
fi
#-------------------------------------------------------------------------------
#  initial checks
#-------------------------------------------------------------------------------
# is mpirun command $MPIRUN set?
if [[ ! -x $MPIRUN ]] ; then 
  echo "error in $0: mpirun command unknown or not executable"
  echo "value of MPIRUN: $MPIRUN"
  exit 1
fi
# is parallel paw executable $PARCODE set?
if [[ ! -x $PARCODE ]] ; then 
  echo "error in $0: parallel PAW executable unknown or not executable"
  echo "value of PARCODE: $PARCODE"
  exit 1
fi
#
###########################################################################
# strip extension .cntl
###########################################################################
JOB=${JOB%.cntl}
#
###########################################################################
# submit job: This may be specific to the local computer
#   use variables
#     - $JOB    rootname of the project
#     - $NNODES
#     - $PARODE name of the ppaw executable, e.g. ppaw_fast.x
#     - $MPIRUN submit command
###########################################################################
export OMP_NUM_THREADS_OLD=$OMP_NUM_THREADS # store previous value for reset
export OMP_NUM_THREADS=1
execute "export OMP_NUM_THREADS=$OMP_NUM_THREADS"
#
# a new TMPDIR is constructed to avoid an error within mpirun due a too long 
#       file name. The directory is deleted and the variable is set  to the
#       original value after execution
export OLDTMPDIR=${TMPDIR}
export TMPDIR=$(mktemp -d /tmp/doppaw.XXXXXXX)
#
execute "${MPIRUN}  -np $NNODES --oversubscribe \
                      $PARCODE \
                      $JOB.cntl 1>out 2>&1 "
export OMP_NUM_THREADS=$OMP_NUM_THREADS_OLD
rmdir ${TMPDIR}
TMPDIR=${OLDTMPDIR}
exit 0

