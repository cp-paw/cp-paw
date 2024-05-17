#!/bin/bash 
###############################################################################
#
#          NAME: paw_scan.sh
#
#         USAGE: paw_scan(.sh) options SRCROOT
#
#       OPTIONS: mandatory: wr optional: hnm0bv
#
#   DESCRIPTION: perform a sequence of paw calculations with different
#                input parameters.
#
#       OPTIONS:
#          -h    print this help message
#          -r "name value1 value2 ...": replacement rule.
#                Every occurance of @name@ will be replaced by one of the 
#                specified values in thefiles with rootname $SRCROOT.
#                Each calculation will be performed in its own directory 
#                $SRCROOT_$value.
#          -w whatid:
#                whatid=run:  Perform CP-PAW simulation. 
#                             Requires .cntl and .strc files
#                whatid=fast: Perform CP-PAW simulation.
#                             Requires .cntl and .strc files
#                whatid=dos:  Make density-of-states. Requires .dcntl file.
#        -n nnodes: submit parallel job with the specified number of processes
#        -j njobs: nr. of paw jobs allowed to run simultaneously (only scalar)
#        -m hostname: limit execution to the specified host
#        -0: dry-run
#        -b: directory containing the paw executables
#        -v: verbose
#
# SRCROOT is the rootname of the source files and the common part of
# the root name for all projects. It must contain the a relative or
# absolute path (i.e. one "/").. The root names of individual projects
# will be formed as TARGETROOT=$SRCROOT_$VALUE/$PROJECT_$VALUE, where
# PROJECT is $SRCROOT stripped of the directory part and VALUE is
# assigned by the rule specified by option -r.
#
# If the target directory does not contain a restart file, a source
# restart file $SRCROOT.rstrt, if present and needed, will be copied
# into the target directory.
#
#  Example:
#       paw_scan.sh -w run -r "EPW 20 30 40" ./si2
#
#   REQUIREMENTS: doppaw.sh, paw_resolve.sh, paw_waittillempty.sh
#                  paw_dos.x
#
#         AUTHOR: Peter E. Bloechl; peter.bloechl@tu-clausthal.de
#
#        CREATED: Dec. 15, 2013
#
###############################################################################
#-------------------------------------------------------------------------------
# help message
#-------------------------------------------------------------------------------
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \tpaw_scan(.sh) options SRCROOT\n"
USAGE="$USAGE \n"
USAGE="$USAGE Purpose:\n"
USAGE="$USAGE \t perform a sequence of paw calculations with different \
                 input parameters."
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -h \t print this help message \n"
USAGE="$USAGE \t -r \"name value1 value2 ...\": replacement rule. \n"
USAGE="$USAGE \t\t  every occurance of @name@ will be replaced by \
                    one of the specified values\n"
USAGE="$USAGE \t\t  in the files with rootname \$SRCROOT. \n"
USAGE="$USAGE \t\t  Each calculation will be performed in its own directory \
                    \$SRCROOT_\$value. \n"
USAGE="$USAGE \t -w whatid: \n"
USAGE="$USAGE \t\t whatid=run perform CP-PAW simulation. \
                   Requires .cntl and .strc files\n"
USAGE="$USAGE \t\t whatid=fast perform CP-PAW simulation\
                   Requires .cntl and .strc files\n"
USAGE="$USAGE \t\t whatid=dos make density-of-states.
                   Requires .dcntl file.\n"
USAGE="$USAGE \t -n nnodes: submit parallel job with the specified number \
                             of processes\n"
USAGE="$USAGE \t -j njobs: nr. of paw jobs allowed to run simultaneously \
                           (only scalar)\n"
USAGE="$USAGE \t -m hostname: limit execution to the specified host\n"
USAGE="$USAGE \t -0: dry-run (creates files but does not run jobs)\n"
USAGE="$USAGE \t -b: directory containing the paw executables\n"
USAGE="$USAGE \t -v: verbose\n"
USAGE="$USAGE \n"
USAGE="$USAGE the argument SRCROOT is the rootname of the source files and\n"
USAGE="$USAGE the common part of the root name for all projects.\n"
USAGE="$USAGE It must contain the a relative or absolute path \
              (i.e. one \"/\")..\n"
USAGE="$USAGE The root names of individual projects will be formed as\n "
USAGE="$USAGE TARGETROOT=\$SRCROOT_\$VALUE/\$PROJECT_\$VALUE.\n"
USAGE="$USAGE where PROJECT is \$SRCROOT stripped of the directory part\n"
USAGE="$USAGE and VALUE is assigned by the rule specified by option -r.\n"
USAGE="$USAGE \n"
USAGE="$USAGE If the target directory does not contain a restart file,\n"
USAGE="$USAGE a source restart file \$SRCROOT.rstrt, if present and needed, \n"
USAGE="$USAGE will be copied into the target directory.\n"
USAGE="$USAGE \n"
USAGE="$USAGE Example:\n"
USAGE="$USAGE \t paw_scan.sh -w run -r \"EPW 20 30 40\" ./si2 \n"
USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
#  individual data
#-------------------------------------------------------------------------------
# name of the bin directory holding the executable ppaw_fast.x or paw_fast.x
# name do not allow any trailing blanks in PAWXDIR!
THISDIR=$(pwd)
export NJOBS=3  # number of paw jobs allowed to run simultaneously (only scalar)
export NNODES=0       # number of nodes / scalar with NNODES=0
export NCOREPERNODE=2  # number of cores per node (2)
export PAWXDIR=$(which paw_fast.x); PAWXDIR=${PAWXDIR%paw_fast.x}
export SELECTEDHOST=""
export DRYRUN=no
export VERBOSE=no
#-------------------------------------------------------------------------------
#  implement dry-run
#-------------------------------------------------------------------------------
function execute(){
    if [[ "$DRYRUN" != yes && "$DRYRUN" != 'no' ]] ; then
      echo "error in $0: illegal value '$DRYRUN' of DRYRUN"
      echo "DRYRUN must be either 'yes' or 'no'"
      exit 1
    fi
    if [[ "$VERBOSE" != yes && "$VERBOSE" != 'no' ]] ; then
      echo "error in $0: illegal value '$VERBOSE' of VERBOSE"
      echo "VERBOSE must be either 'yes' or 'no'"
      exit 1
    fi
    # || is "or" in [[...]], && is "and" in [[..]] 
    if [[ "${DRYRUN}" = "yes" || "${VERBOSE}" = "yes" ]] ; then
      echo "${@}"
    fi
    if [[ "${DRYRUN}" = "no" ]] ; then
      eval "$@"
      export PID=$!
    fi
}
#-------------------------------------------------------------------------------
#  count the number of running jobs with PID in PIDS
#-------------------------------------------------------------------------------
function countjobs {
################################################################################
## returns the number of active jobs on a list PIDS on project id-s           ##
###############################################################################
  local NACTIVE=0
  local X
  for X in ${PIDS} ; do
    if [[ -n  $(ps -p $X -o pid= ) ]] ; then  NACTIVE=$(( $NACTIVE + 1 )); fi
  done
echo $NACTIVE
}
#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
while getopts :hv0w:r:b:x:n:m:f:j: OPT  ; do
  echo OPT=$OPT
  case $OPT in
    r)  #rule:   "id val1 val2 val3 ....". Use apostrophes!
      RULE="$OPTARG"
      PAR1NAME=""
      PAR1VALS=""
      for X in $RULE ; do
        if [[ -z $PAR1NAME ]] ; then
          PAR1NAME=${X}
        else
          PAR1VALS="$PAR1VALS $X"
        fi
      done
      ;;
    w)   # what: run|fast|dos) 
      WHAT=$OPTARG
      ;;
    n)   #nodes:
      NNODES=$OPTARG
      ;;
    j)   #jobs
      NJOBS=$OPTARG
      ;;
    m)   #machine:
      SELECTEDHOST=$OPTARG
      ;;
    b)   #builddirectory:
      PAWXDIR=$OPTARG
      ;;
    0)   #nothing:
      DRYRUN=yes
#      set -n
      ;;
    v)   #verbose
      VERBOSE=yes
#      set -v
#      set -x
      ;;
    h)   # help
      echo -e $USAGE
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
SRCROOT=$1  #  Directory containing all subprojects relative to THISDIR
#
#===================================================================
# check arguments
#===================================================================
if [[ -z $SRCROOT ]] ; then 
  echo "error in $0: missing argument SRCROOT" >&2
  exit
fi
SRCNAME=${SRCROOT##*/}
SRCDIR=${SRCROOT%/*}
#
#
#
if [[ -z $PAR1NAME ]] ; then 
  echo  "error in $0: missing or incomplete argument -r" >&2
  echo  "rule=$RULE" >&2
  echo  "PAR1NAME=$PAR1NAME" >&2
  echo  "PAR1VALS=$PAR1VALS" >&2
  exit
fi

if [[ -z $WHAT ]] ; then 
  echo "error in $0: option -w not specified" >&2
  exit 1
fi
#
#
#
case $WHAT in
  fast|run)
    EXECTBLE=paw_fast.x
    if [[ $NNODES != 0 ]] ; then
      EXECTBLE=ppaw_fast.x
    fi
    CNTL=cntl
    FILES="cntl strc"
    # needs strc 
    WHAT=run
    ;;
  dos)
    EXECTBLE=paw_dos.x
    CNTL=dcntl
    FILES="dcntl"
    # needs strc_out? pdos 
   ;;
  *) 
    echo "error in $0: illegal value $WHAT for WHAT" >&1
    exit 1
    ;;
esac
EXECTBLE=$PAWXDIR/$EXECTBLE
#
# run only on selected host
#
if [[ -n $SELECTEDHOST ]] ; then 
  if [[ ${HOSTNAME%%.*} != $SELECTEDHOST ]] ; then
    if [[ $DRYRUN = no ]] ; then
      echo "error in $0: not on $SELECTEDHOST" >&2
      echo "actual host name: ${HOSTNAME%%.*}" >&2
      exit 1
    fi
  fi
fi
#
#  variables defined: 
#    SRCROOT
#    SELECTEDHOST EXECTBLE CNTL
#    PAR1NAME PAR1RULES
#    SRCROOT
#    DRYRUN VERBOSE
#
#-------------------------------------------------------------------------------
#    report arguments
#-------------------------------------------------------------------------------
if [[ $VERBOSE = yes ]] ; then
  echo ========================================================
  echo "WHAT........: $WHAT"
  echo "SRCROOT.....: $SRCROOT"
  echo "SRCNAME.....: $SRCNAME"
  echo "SRCDIR......: $SRCDIR"
  echo "FILES.......: $FILES"
  echo "PAR1NAME....: $PAR1NAME"
  echo "PAR1VALS....: $PAR1VALS"
  echo "executable..: $EXECTBLE"
  echo "control file: $CNTL"
  echo "HOSTNAME....: $HOSTNAME"
  echo "Nr. Jobs....: $NJOBS"
  echo "========================================================"
fi
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    initialize: make directories and copy files
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
export PIDS="" # list of project id's
#
for PAR1VAL in $PAR1VALS; do
  TARGETNAME=${SRCNAME}'_'${PAR1VAL}
  TARGETDIR=${SRCDIR}/${TARGETNAME}
  TARGETROOT=${TARGETDIR}/$TARGETNAME
  if [[ $VERBOSE = yes ]] ; then
    echo "========================================================"
    echo "parval.......: ${PAR1VAL}"
    echo "targetname...: ${TARGETNAME}"
    echo "targetroot...: ${TARGETROOT}"
    echo "targetdir....: ${TARGETDIR}"
  fi
  #
  # create target directory
  #
  if [[ ! -d $TARGETDIR ]] ; then mkdir $TARGETDIR; fi
  #
  #
  for EXT in $FILES; do
     SOURCE=${SRCROOT}.$EXT
     TARGET=${TARGETROOT}.$EXT
     if [[ ! -e ${SOURCE} ]] ; then 
       echo "error in $0: missing file ${SOURCE}" >&2; exit 1
     fi
     #
     # remove target if older than source
     #
     if [[ -e $TARGET ]] ; then
       if [ ${TARGET} -ot ${SOURCE} ] ; then rm $TARGET; fi
     fi
     #
     # copy if target does not exist or is older than source
     #
     if [[ ! -e ${TARGET} ]] ; then
       echo copying file into $TARGET
       paw_resolve.sh -r "${PAR1NAME}=$PAR1VAL" \
                      -r "SRCROOT=$SRCROOT" \
                      -i "${SOURCE}" \
                      -o "${TARGET}"
       RC=$?
       if [[ $RC -ne 0 ]] ; then
         echo "error in $0: paw_resolve.sh failed to expand ${SOURCE}" >&2
         exit 1
       fi
     fi
  done
  #
  #   copy restart file
  #
  if [[ $WHAT = run ]] ; then
    SOURCE=${SRCROOT}.rstrt  
    TARGET=${TARGETROOT}.rstrt  
    if [[ ! -e ${TARGET} ]] ; then  # copy only if target does not exist
      if [[ -e ${SOURCE} ]] ; then # copy only if source exists
        echo copying file into ${TARGET}
        execute "cp ${SOURCE} ${TARGET}"
        echo 'restart file copied...' `date`
      fi
    fi
  fi
  #
  # EXECUTE
  #
  cd $TARGETDIR
  if [[ $NNODES = 0 ]] ; then
    # pause if the number of active jobs exceeds maximum
    while [[ $(countjobs) -ge ${NJOBS} ]] ; do sleep 2; done
    #
#    execute "paw_waittillempty.sh -n $NJOBS"
    execute "$EXECTBLE ${TARGETNAME}.$CNTL 1>out 2>&1 &"
    PIDS="$PIDS $PID"   # store the project id returned by 'execute'
  else
    if [[ WHAT != run ]] ; then
      echo "error in $0: only -w run can be executed parallel" >&2
      echo "WHAT (-w)=$WHAT" >&2
      echo "number of nodes requested: $NNODES" >&2
      exit 1
    fi
    execute "doppaw.sh ${TARGETNAME} ${NNODES}"
    PIDS="$PIDS $PID"   # store the project id returned by 'execute'
  fi
  cd $THISDIR
done
#
#============================================================================
#==  wait until all submitted jobs ($PIDS) are finished
#============================================================================
while [[ $(countjobs) -gt 0 ]] ; do 
  sleep 10 
  echo "time=$(date "+ %H:%M:%S"); waiting for $(countjobs) to finish...."
done
#
exit 0
