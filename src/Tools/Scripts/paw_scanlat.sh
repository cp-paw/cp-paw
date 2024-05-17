#!/bin/bash
###############################################################################
#
#          NAME: paw_scanlat.sh
#
#         USAGE: paw_scanlat(.sh) options 
#
#       OPTIONS: mandatory: p optional: hlubej0v 
#
#   DESCRIPTION: perform a sequence of paw calculations with different
#                lattice constants or updates data files with its results.
#                produces the file latscan.dat, which lists the resuling
#                volumes, total energies and percent changes of the 
#                lattice constant
#
#       OPTIONS:
#          -h    print this help message
#          -p    rootname of the project to be considered (mandatory)
#          -l    list of percent changes of the lattice constant
#                 (enter in apostrophes such as "98 100 102")
#                default: "96 97 98 99 100 101 102 103 104"
#          -u    toggle: if -u is present, only update of data files, but no paw simulations
#                if -u is absent the calculations are performed, but not updated
#          -e    paw command (default paw_fast.x)
#          -j    nr. of paw jobs allowed to run simultaneously (default 8)
#          -0    dry-run (reconstruct *.dat files only)
#          -v    verbose
#
#  prerequisites: nio.cntl 
#                 nio.strc 
#        (nio is to be replaced by the name of the substance )
#
#  performs set of calculations for different lattice constants.
#
#  1) starts from a rootname ROOT specified as project by the option
#     -p.  ROOT may contain the full path or the relative path with
#     respect to the current directory of the parent process. As
#     precondition there is a file $ROOT.cntl and a file $ROOT.strc
#
#  2) with the option -l a list of percentage changes of the lattice
#     constant is provided
#
#  3) For each percentage change $X of the lattice constant in the
#     list, a calculation is performed in its own directory
#     ${DIROFX}=${ROOT}_$X. The rootname of this project is
#     ${ROOTOFX}=${DIROFX}/${DIROFX##*/}.
#
#  nio.strc must contain exactly one line containing "LUNIT=" 
#  followed by the lattice constant in atomic units, or "LUNIT[AA]="
#  followed by the lattice constant in angstrom.
#  This value of LUNIT will be replaced by the modified values for the 
#  individual runs.
#
#  specify k-points in the structure file with !STRUCTURE!KPOINTS:DIV=
#  rather than using R=. Otherwise the k-point grid changes between 
#  different values of lattice constant and the restart file becomes 
#  unusable.
#
#  result       etot.dat 
#               gap.dat
#               murn.in
#
#  Example:
#       paw_scanlat.sh -p nio -j3 -l "96 98 100 102 104" ./si2
#
#   REQUIREMENTS: paw_get.sh, paw_fast.x
#
#         AUTHOR: Peter E. Bloechl; peter.bloechl@tu-clausthal.de
#
#        CREATED: Sept. 2, 2014
#
################################################################################
#==============================================================================
# initialize variables and scan argument list
#==============================================================================
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \t paw_scanlat[.sh] options\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -p projectname of the project to be considered (mandatory)\n"
# assumption: project name does not contain the directory
USAGE="$USAGE \t -l list of percent changes of the lattice constant\n"
USAGE="$USAGE \t    (enter in apostrophes \"...\")\n"
USAGE="$USAGE \t     \"96 97 98 99 100 101 102 103 104\"\n"
USAGE="$USAGE \t -u updates result data files (no paw simulations)\n"
USAGE="$USAGE \t -e paw command (default: paw_waittillempty.sh -n 8;\
                                                              paw_fast.x)\n"
USAGE="$USAGE \t -j nr. of paw jobs allowed to run simultaneously\n"
USAGE="$USAGE \t -0 dry-run (reconstruct *.dat files only)\n"
USAGE="$USAGE \t -v verbose\n"
USAGE="$USAGE \t -h print this help message \n"
USAGE="$USAGE \n"
USAGE="$USAGE preconditions: \n" 
USAGE="$USAGE \t projectname.cntl \n" 
USAGE="$USAGE \t projectname.strc \n" 
USAGE="$USAGE \n"
#
export THISDIR=$(pwd)   # current directory
# PAWXDIR is the directory holding the paw executables
export PAWXDIR=$(which paw_fast.x); PAWXDIR=${PAWXDIR%paw_fast.x}
#-------------------------------------------------------------------------------
#  implement dry-run
#-------------------------------------------------------------------------------
function execute(){
    if [[ "$DRYRUN" != yes && "$DRYRUN" != 'no' ]] ; then
      echo "error in $0" >&2
      echo "illegal value '$DRYRUN' of DRYRUN" >&2
      echo "DRYRUN must be either 'yes' or 'no'" >&2
      exit 1
    fi
    if [[ "$VERBOSE" != yes && "$VERBOSE" != 'no' ]] ; then
      echo "error in $0" >&2
      echo "illegal value '$VERBOSE' of VERBOSE" >&2
      echo "VERBOSE must be either 'yes' or 'no'" >&2
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
#  strip leading and trailing whitespace from a bash variable
#-------------------------------------------------------------------------------
function stripwhitespace {
##########################################################################
# strips leading and trailing space characters from the variable
# passed as argument. 
#
# -- [:space:] is a posix character class which includes tab, newline,
#  vertical tab, form feed, carriage return, and space.
# -- *([:space:]) is a regular expression for zero or more spaces.
# -- ${X#RegEx} removes the shortest leading match of the regular
#  expression RegEx in the variable X
# -- ${X%RegEx} removes the shortest trailing match of the regular
#  expression RegEx in the variable X
##########################################################################
local X=$1    # use all arguments as one string. 
X=${X##*([:space:])} # remove leading whitespace characters
X=${X%%*([:space:])} # remove trailing whitespace characters
echo "$X"
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
#
#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
export NAME # root name of the project including relative of absolute path
ALATLIST="96 97 98 99 100 101 102 103 104"
DRYRUN=no
VERBOSE=no
UPDATE=no
NJOBS=8
while getopts :h0e:p:l:j:uv OPT ; do
  case $OPT in
    p)
      ROOT=$OPTARG
      ;;
    l)
      ALATLIST=$OPTARG
      ;; 
    e)
      PAWCMD=$OPTARG
      ;;
    u)
      UPDATE=yes
      ;;
    j)
      NJOBS=$OPTARG
      ;;
    0)
      DRYRUN=yes
      ;;
    v)
      VERBOSE=yes
      ;;
    h)
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
if [[ -e $1 ]] ; then echo "error in $0: argument present. none expected" >&2; fi
if [[ -z $ROOT ]] ; then
  echo "error in $0: mandatory option -p missing" >&2
  exit 1
fi
if [[ ! -f $ROOT.strc ]] ; then
  echo "error in $0: file $ROOT.strc is missing" >&2
  exit 1
fi
if [[ ! -f $ROOT.cntl ]] ; then
  echo "error in $0: file $ROOT.cntl is missing" >&2
  exit 1
fi
#
# when PAWCMD has not been defined via arguments, set default
if [[ -z $PAWCMD ]] ; then
  PAWCMD="paw_fast.x"
fi
#
#   report setting of variables
#   
if [ $VERBOSE = yes ] ; then
  echo "root name of the project...................: ${ROOT}"
  echo "list of percent changes of lattice constant: ${ALATLIST}"
  echo "paw command................................: ${PAWCMD}"
  echo "dry-run....................................: ${DRYRUN}"
  echo "update.....................................: ${UPDATE}"
  echo "max. nr. of jobs running simultaneously....: ${NJOBS}"
fi
#
#==============================================================================
# Divide root into a directory opart and a non-directory part
#==============================================================================
export NAME=${ROOT##*/}     # non-directort  part of $ROOT
export DIR=${ROOT%${NAME}}  # directory part of $ROOT
if [[ -z ${DIR} ]] ; then 
  DIR=${THISDIR} 
else
  DIR=${THISDIR}/${DIR##${THISDIR}}
fi
ROOT=${DIR}/${NAME}
if [ $VERBOSE = yes ] ; then echo DIR =\"${DIR}\" ; fi
if [ $VERBOSE = yes ] ; then echo NAME=\"${NAME}\" ; fi
if [ $VERBOSE = yes ] ; then echo ROOT=\"${ROOT}\" ; fi
#
#==============================================================================
# pick out the reference lattice constant from the structure file
#==============================================================================
LINE=$(grep -i 'LUNIT' $ROOT.strc)
echo "LINE(1) ${LINE}"
LINE=${LINE#*[Ll][Uu][Nn][Ii][Tt]} #strip everything up to LUNIT
echo "LINE(2) $LINE"
if [[ -n ${LINE%%\[AA\]=*} ]] ; then
  echo unit is bohr
  AA=$(echo "1. / 0.529177" | bc -l )   # conversion factor into atomic units
else
  echo unit is angstrom
  AA=$(echo "1." | bc -l )
fi
#echo "LINE(3) $LINE"
LUNIT=${LINE#*=} # remove all till next preceeding equal sign (such as [AA])
#echo "LUNIT(1) --$LUNIT--"
LUNIT=$(stripwhitespace "$LUNIT") # stripwhitespace is a function
                                  # defined above. Keep the
                                  # apostrophes around $LUNIT, so that
                                  # it is interpreted as one string
                                  # with spaces rather than an array
                                  # of strings separated by spaces
LUNIT=${LUNIT%%[[:space:]]*} # remove trailing text
ALAT0=$LUNIT
if [ $VERBOSE = yes ] ; then echo ALAT0=\"${ALAT0}\"; fi
#
#==============================================================================
# loop crystal calculation through different lattic constants                ==
#==============================================================================
if [[ $UPDATE = no ]] ; then
  export PIDS="" # list of project id's
  for X in $ALATLIST; do
    if [[ $VERBOSE = yes ]] ; then echo doing $X ; fi 
    export DIROFX=${DIR}/${NAME}_$X
    export ROOTOFX=${DIROFX}/${NAME}_$X
    export ALATX=$(echo "$X /100 * $ALAT0 " | bc -l)
    if [[ $VERBOSE = yes ]] ; then echo DIROFX=${DIROFX} ; fi 
    if [[ $VERBOSE = yes ]] ; then echo ROOTOFX=${ROOTOFX} ; fi 
    if [[ $VERBOSE = yes ]] ; then echo ALATOFX=${ALATX} ; fi 
    if [[ ! -d $DIROFX ]] ; then mkdir ${DIROFX} ; fi

    #===========================================================================
    # construct structure control file from template by replacing alat
    #===========================================================================
    execute "sed \"/[Ll][Uu][Nn][Ii][Tt]/s/${ALAT0}/${ALATX} /g\" \
                                         ${ROOT}.strc >${ROOTOFX}.strc"
    #==========================================================================
    # run paw_fast.x and place data "gap.dat" and "etot.dat"
    #==========================================================================
    execute "cp ${ROOT}.cntl ${ROOTOFX}.cntl"
    if [[ ! -e ${ROOTOFX}.rstrt ]] ; then
      if [[ -e ${ROOT}.rstrt ]] ; then
        execute "cp ${ROOT}.rstrt ${ROOTOFX}.rstrt"
      else
        # test, whether restart file is required
        if [[ -z $(grep -i START=*T  ${ROOTOFX}.cntl) ]]; then 
          echo "source restart file ${ROOTOFX}.rstrt does not exist"
          exit 1
        fi
      fi
    fi
    #
    # pause if the number of active jobs exceeds maximum
    while [[ $(countjobs) -ge ${NJOBS} ]] ; do sleep 2; done
    #
    execute "cd ${DIROFX}"
    execute "${PAWCMD} ${ROOTOFX}.cntl 1>${ROOTOFX}.out 2>&1 &" 
    PIDS="$PIDS $PID"   # store the project id returned by 'execute'
    execute "cd ${THISDIR}"
  done
fi
#
#============================================================================
#==  wait until all submitted jobs ($PIDS) are finished
#============================================================================
while [[ $(countjobs) -gt 0 ]] ; do 
  sleep 10 
  echo "time=$(date "+ %H:%M:%S"); waiting for $(countjobs) to finish...."
done

# export NACTIVE=1
# while [ ${NACTIVE} -gt 0 ] ; do
#   NACTIVE=0
#   for X in ${PIDS} ; do
#     if [[ -n  $(ps -p $X -o pid= ) ]] ; then  NACTIVE=$(( $NACTIVE + 1 )); fi
#   done
#   echo "No.Jobs=$NACTIVE; time=$(date "+ %H:%M:%S"); waiting...."
#   sleep 2
# done
#
#==============================================================================
# update data files
#==============================================================================
if [[ $UPDATE = yes ]] ; then
  if [[ "$DRYRUN" = "no" ]] ; then
    echo "# VOLUME[cube-abohr] ETOT[H] ${X}[percent]" >${DIR}/latscan.dat
  fi
  for X in $ALATLIST; do
    if [[ $VERBOSE = yes ]] ; then echo collecting $X ; fi 
    export DIROFX=${DIR}/${NAME}_$X
    export ROOTOFX=${DIROFX}/${NAME}_$X
    export ALATX=$(echo "$X /100 * $ALAT0 " | bc -l)
    if [[ $VERBOSE = yes ]] ; then echo DIROFX=${DIROFX} ; fi 
    if [[ $VERBOSE = yes ]] ; then echo ROOTOFX=${ROOTOFX} ; fi 
    if [[ $VERBOSE = yes ]] ; then echo ALATOFX=${ALATX} ; fi 
    # "paw_get.sh -w volume" takes the last lattice vectors from the prot file
    VOLUME=$(paw_get.sh -n -w volume -u a0^3 ${ROOTOFX})
    ETOT=$(paw_get.sh -n -w etot -u H ${ROOTOFX})
    if [[ $VERBOSE = yes ]] ; then echo VOLUME=${VOLUME} ; fi 
    if [[ $VERBOSE = yes ]] ; then echo ETOT=${ETOT} ; fi 
    if [[ "$DRYRUN" = "no" ]] ; then
      echo ${VOLUME} ${ETOT} ${X} >>${DIR}/latscan.dat
    fi
    #
    #  clean up
    # 
    if [[ "$DRYRUN" = "no" ]] ; then
      LIST="${ROOTOFX}_constr.report ${ROOTOFX}_r.tra ${ROOTOFX}.pdos \
            ${ROOTOFX}.strc_out "
      for Y in $LIST ; do
        if [[ $VERBOSE = yes ]] ; then echo cleaning $Y ; fi
        if [[ -e $Y ]] ; then rm $Y ; fi
      done
      for Y in ${ROOTOFX}_stpfor*.myxml; do
        if [[ -e $Y ]] ; then rm $Y ; fi
      done
    fi
  done
fi
execute "cd ${THISDIR}"
exit 0
