#!/bin/bash 
###############################################################################
#
#          NAME: paw_do.sh
#
#         USAGE: paw_do(.sh) options 
#
#       OPTIONS: optional: hv0fsdw 
#
#   DESCRIPTION: frontend for paw calculations or tools.
#
#       OPTIONS:
#          -h    print this help message
#          -0:   dry-run
#          -v:   verbose
#          -R root  specify root file. (Default is to find it.)
#          -F    execute paw_fast.x
#          -K    soft kill of simulation
#          -P    execute doppaw.sh (ppaw_fast.x)
#          -D    execute paw_dos.x and paw_dosplot.x
#          -E arg interface for use in shell scripts.   
#          -B    execute paw_bands.x and xmgrace -nxy bands.dat
#          -S    execute paw_strc.x
#          -W    execute paw_wave.x
#          -T    execute paw_tra.x
#          -A    execute paw_show.sh -ce
#          -Z    execute tail -f ROOT.prot
#       uppercase options are chosen to select a particular tool. Lowercase
#       options (except hv0) are treated as options for the preceeding tool 
#
#  Example:
#       paw_do.sh -F 
#
#   REQUIREMENTS: 
#                  doppaw.sh
#
#                  paw_dos.x,paw_dosplot.x,paw_show.sh,paw_fast.x,paw_strc.x
#                  paw_wave.x,paw_tra.x,paw_bands.x
#
#          the scripts rely on the fact that exactly one file with the
#          required extension lies in the local directory.
#
#         AUTHOR: Peter E. Bloechl; peter.bloechl@tu-clausthal.de
#
#        CREATED: Nov. 17, 2014
#
###############################################################################
#-------------------------------------------------------------------------------
# help message
#-------------------------------------------------------------------------------
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \tpaw_do options\n"
USAGE="$USAGE \n"
USAGE="$USAGE Purpose:\n"
USAGE="$USAGE \t frontend for paw code and tools"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -h \t print this help message \n"
USAGE="$USAGE \t -A: paw_show.sh -ce ROOT\n"
USAGE="$USAGE \t -B: paw_bands.x ROOT.bcntl\n"
USAGE="$USAGE \t\t   xmgrace -nxy bands.dat\n"
USAGE="$USAGE \t -F: paw_fast.x ROOT.cntl 1>out 2>&1 & \n"
USAGE="$USAGE \t -D: paw_dos.x ROOT.dcntl \n"
USAGE="$USAGE \t\t   paw_dosplot.x ROOT.dpcntl; xmgrace -batch ROOT.bat\n"
USAGE="$USAGE \t -K: touch ROOT.exit (regular stop of simulation) \n"
USAGE="$USAGE \t -P: doppaw.sh -n 6 ROOT (parallel ppaw_fast.x with 6 tasks)\n"
USAGE="$USAGE \t -R  root:  specify root file. (Default is to find it.)\n"
USAGE="$USAGE \t -S: paw_strc.x -c ROOT\n"
USAGE="$USAGE \t -T: paw_tra.x ROOT.tcntl\n"
USAGE="$USAGE \t -W: paw_wave.x ROOT.wcntl\n"
USAGE="$USAGE \t -Z: tail -f ROOT.prot\n"
USAGE="$USAGE \t -0: dry-run (creates files but does not run jobs)\n"
USAGE="$USAGE \t -v: verbose\n"
USAGE="$USAGE \n"
USAGE="$USAGE Example:\n"
USAGE="$USAGE \t paw_do -F \n"
USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
#  individual data
#-------------------------------------------------------------------------------
# name of the bin directory holding the executable ppaw_fast.x or paw_fast.x
# name do not allow any trailing blanks in PAWXDIR!
THISDIR=$(pwd)
export PAWXDIR=$(which paw_fast.x); PAWXDIR=${PAWXDIR%paw_fast.x}
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
    fi
}
#-------------------------------------------------------------------------------
#  findroot
#-------------------------------------------------------------------------------
# searches the current directory for a file with a specified extension
# and returns the file name if it is unique. Otherwise it returns an empty 
# string. 
# use as FILE=$(findroot cntl)
function findroot(){
    local LIST=$(echo *.$1)
    local RESULT=
    local X
    if [[ -n ${ROOT} ]] ; then
      if [[ -f ${ROOT}.$1 ]] ; then
        RESULT=${ROOT}.$1  # use variable root if defined, and if file exists
      else
        RESULT=
      fi
    else
      for X in $LIST; do
        if [[ -z $RESULT ]] ; then
          RESULT=$X
        else  # return nothing, if more than one match!
          RESULT=
          LIST=
        fi
      done
    fi
    echo $RESULT
}

#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
export ROOT=
#
OPTSTRING=":hv0ABE:FKPR:SDWTZ"
OPTIND=0
while getopts "${OPTSTRING}" OPT  ; do
  case $OPT in
    R) ROOT="$OPTARG" #set explicit root name
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
if [[ "${VERBOSE}" == "yes" ]] ; then
  echo verbose=$VERBOSE
  echo dryrun=$DRYRUN
  echo OPTSTRING=$OPTSTRING
  echo args="$1"
fi
#
#
#
#
OPTIND=0  # reset argument position
while getopts "${OPTSTRING}" OPT  ; do
 case $OPT in
    A) # paw_show.sh -ce $ROOT
      Y=$(findroot prot)
      Y=${Y%.prot}
      if [[ -n $Y ]] ; then      
        execute "paw_show.sh -ce $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.prot file not unique or not present" >&2
        exit 1
      fi
      ;;
    B) # paw_bands.x $ROOT.bcntl
      Y=$(findroot bcntl)
      if [[ -n $Y ]] ; then      
        execute "paw_bands.x  $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.bcntl file not unique or not present" >&2
        exit 1
      fi
      if [[ -e bands.dat ]] ; then      
        execute "xmgrace -free -noask -nxy bands.dat"
      else
        echo "error in $0 for option $OPT" >&2
        echo "file bands.dat not present" >&2
        echo "in the bcntl file, direct output into bands.dat" >&2
        echo "or inspect the file specified in .bcntl with xmgrace." >&2
        exit 1
      fi
      ;;
    D) # paw_dos.x $ROOT.dcntl
      Y=$(findroot dcntl)
      if [[ -n $Y ]] ; then      
        execute "paw_dos.x  $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.dcntl file not unique or not present" >&2
        exit 1
      fi
      Y=$(findroot dpcntl)
      if [[ -n $Y ]] ; then      
        execute "paw_dosplot.x $Y"
        Y=${Y%.dpcntl}
        execute "xmgrace -free -noask -batch $Y.bat"
      fi
      ;;
    E) # execute:  entry for use in shell scripts
      Y=$(findroot cntl)
      Y=${Y%.cntl}
      if [[ -z $Y ]] ; then      
        echo "error in $0 for option $OPT" >&2
        echo "*.cntl file not unique or not present" >&2
        exit 1
      fi
      case $OPTARG in
        P)
          execute "nohup doppaw.sh -n 6 $Y"
          ;;
        F)
          execute "paw_fast.x  $Y 1>out 2>&1"
          ;;
        *)
          echo "error in $0 for option $OPT" >&2
          echo "option argument $OPTARG not recognized" >&2
          exit 1
          ;;
      esac
      ;;
    F) # paw_fast.x $ROOT.cntl
      Y=$(findroot cntl)
      if [[ -n $Y ]] ; then      
        execute "paw_fast.x  $Y 1>out 2>&1 &"
        execute "while [ ! -e ${Y%.cntl}.prot ] ; do sleep 1; done"
        execute "tail -f ${Y%.cntl}.prot"
       else
        echo "error in $0 for option $OPT" >&2
        echo "*.cntl file not unique or not present" >&2
        exit 1
      fi
      ;;
    K) # touch $ROOT.exit
      Y=$(findroot cntl)
      if [[ -n $Y ]] ; then      
        execute "touch ${Y%.cntl}.exit"
       else
        echo "error in $0 for option $OPT" >&2
        echo "*.cntl file not unique or not present" >&2
        exit 1
      fi
      ;;
    P) # ppaw_fast.x $ROOT.cntl 
      Y=$(findroot cntl)
      Y=${Y%.cntl}
      if [[ -n $Y ]] ; then      
        execute "nohup doppaw.sh -n 10 $Y &"
        execute "while [ ! -e ${Y%.cntl}.prot ] ; do sleep 1; done"
        execute "tail -f ${Y%.cntl}.prot"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.cntl file not unique or not present" >&2
        exit 1
      fi
      ;;
    S) # paw_strc.x -c 
      Y=$(findroot strc_out)
      Y=${Y%.strc_out}
      if [[ -n $Y ]] ; then      
        execute "paw_strc.x -c $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "ROOT not unique or no file *.strc_out present" >&2
        exit 1
      fi
      ;;
    W) # paw_wave.x $ROOT.wcntl
      Y=$(findroot wcntl)
      if [[ -n $Y ]] ; then      
        execute "paw_wave.x  $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.wcntl file not unique or not present" >&2
        exit 1
      fi
      ;;
    T) # paw_tra.x $ROOT.tcntl
      Y=$(findroot tcntl)
      if [[ -n $Y ]] ; then      
        execute "paw_tra.x  $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.tcntl file not unique or not present" >&2
        exit 1
      fi
      ;;
    Z) # tail -f $ROOT.prot
      Y=$(findroot prot)
      if [[ -n $Y ]] ; then      
        execute "tail -f $Y"
      else
        echo "error in $0 for option $OPT" >&2
        echo "*.prot file not unique or not present" >&2
        exit 1
      fi
      ;;
  esac
done
exit 0
