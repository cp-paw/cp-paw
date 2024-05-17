#!/bin/bash
#
# waits until there are less than $1 processes with "paw" in its name running
# This is used to delay startin of a paw_job until there is sufficient room 
# (# of cores) available. With an argument the number of jobs is limited to
# the number given as input argument
#
# check which operating system and analyze computer
#
export USAGE="Usage of $0 \n"
export USAGE="$USAGE \n"
export USAGE="$USAGE paw_waittillempty[.sh] options\n "
export USAGE="$USAGE \n"
export USAGE="$USAGE pauses until the number of paw-related jobs \n"
export USAGE="$USAGE falls below njobx\n"
export USAGE="$USAGE \n"
export USAGE="$USAGE njobx=max number of allowed paw jobs\n"
export USAGE="$USAGE \n"
export USAGE="$USAGE options:\n"
export USAGE="$USAGE -h|?: help\n"
export USAGE="$USAGE -n njobx: specify the maximum of allowed paw jobs\n"
export USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
while getopts :n:h? OPT; do
  case $OPT in
    n)
      NJOBX=$OPTARG
      ;;
    h)
      echo -e $USAGE
      exit 1
      ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0" >&2
      echo "invalid option -$OPTARG" >&2
      echo "retrieve argument list with:" >&2
      echo "$0 -h" >&2
      exit 1
      ;;
    :)
      echo "error in $0" >&2
      echo "option -$OPTARG requires an additional argument" >&2
      exit 1
      ;;
  esac
done
echo NJOBX=$NJOBX
#
#  determine the number of cores as default for NJOBX
#
if [ -z $NJOBX ] ; then   # number of argmuments equals 1?
  export OS=$(uname -s)
  if [ $OS = Darwin ] ; then   # Darwin identifies Apple OSX
    export NCORE=$(sysctl hw.ncpu | awk '{print $2}')
  elif [ $OS = Linux ] ; then 
#    export NCORE=$(lscpu | grep 'CPU(s):' | grep -v ' CPU(s):*' \
#                         | awk -F: '{print $2}')
     if [[ -n "$(which lscpu)" ]] ; then
       export NCORE=$(lscpu | grep 'CPU(s):' | head -1 | cut -d ":" -f2)
     else
        export NCORE=$(grep -c processor /proc/cpuinfo)
     fi
  else
    echo "operating system $OS not recognized" >&2
    exit 1
 fi
  export NCOREPERJOB=2    #assume one virtual per real core
  export NJOBX=$(echo $NCORE  / $NCOREPERJOB | bc)
  echo $NCORE cores on the processor
  echo a max of $NCOREPERJOB cores is used per job, thus 
fi
echo no more than $NJOBX jobs shall be run at a time
#
# pause until one more job fits
#
NJOBS=$(ps -e | grep -epaw_[a-z]*.x | wc -l)
echo $NJOBS jobs are currently running. date: $(date) 
echo entering waiting loop .....
sleep  2  # sleep 2 seconds to allow jobs to start
while [ $NJOBS -ge $(echo $NJOBX | bc) ] ; do
  sleep  10
  NJOBS=$(ps -e | grep -epaw_[a-z]*.x | wc -l)
  echo No.Jobs=$NJOBS; time=$(date "+ %H:%M:%S") waiting....
done
echo $NJOBS jobs are currently running. date: `date` 
echo .... waiting loop finished.  
echo waittillempty terminated.

