#!/bin/bash
############################################################################
#
#  paw_get.sh
#
#  command: paw_get.sh options ROOTNAME
#
#  extracts the last value for etot|gap|homo|lumo|efermi from the
#  specified protocol file and prints either text or number
#
############################################################################
#
#  set up help message $USAGE
#
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \t paw_get.sh options rootname\n"
USAGE="$USAGE extracts data from file rootname.prot\n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -w whatid\n"
USAGE="$USAGE \t\t whatid=etot  \t prints total energy in Hartree (default)\n"
USAGE="$USAGE \t\t whatid=gap   \t prints fundamental band gap in eV\n"
USAGE="$USAGE \t\t whatid=efermi\t prints fermi level eV\n"
USAGE="$USAGE \t\t whatid=homo  \t prints energy of the highest occupied \
                                          state in eV\n"
USAGE="$USAGE \t\t whatid=lumo  \t prints energy of the lowest unoccupied \
                                          state in eV\n"
USAGE="$USAGE \t\t whatid=volume \t prints volume of the unit cell \n"
USAGE="$USAGE \t\t whatid=epw    \t prints wavefunction plane wave cutoff \n"
USAGE="$USAGE \t\t whatid=epwrho \t prints density plane wave cutoff \n"
USAGE="$USAGE \t\t whatid=nkpt   \t prints the number of k-points\n"
USAGE="$USAGE \t\t whatid=kdiv   \t prints the divisions for the k-grid\n"
USAGE="$USAGE \t -n prints the result as number or according to -u\n"
USAGE="$USAGE \t -t prints the result as text (default)\n"
USAGE="$USAGE \t -u unit: 1|ev|kj/mol|kcal/mol|h|ry|angstrom^3|aa^3|cubangstrom|a0^3|abohr^3|cubabohr \n"
USAGE="$USAGE \t -v verbose\n"
USAGE="$USAGE \t -h prints this help message\n"
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
#  implement unit conversion
#-------------------------------------------------------------------------------
function convertunit(){
  # takes a string as argument containing a number and a unit string.
  # returns a new variable "AUVALUE" in hartree atomic units
  local NEWUNIT=$1
  local ARGSTRING=$2
  local ARGVALUE=$(echo "$ARGSTRING" | awk '{print $1}')
  local ARGUNIT=$(echo "$ARGSTRING" | awk '{print $2}')
  #
  local HARTREE=1.
  local RYDBERG=0.5
  local EV=$(echo "scale=8; 1. / 27.211 " | bc -l)
  local KJBYMOL=$(echo "scale=8; 1. / 2625.500223430069 " | bc -l)
  local KCALBYMOL=$(echo "scale=8; 1. / 627.5096463920391 " | bc -l)
  local ANGSTROM=$(echo "scale=8; 1. / 0.529177 " | bc -l)
  local CUBANGSTROM=$(echo "scale=8; $ANGSTROM ^ 3 " | bc -l)
  local ABOHR=1.
  local CUBABOHR=1.
  #
  local NEWVALUE
  # check if $ARGVALUE is a number
  if [[ $(echo "${ARGVALUE}" \
         | grep -x -e "[-]\?[0-9]*[.]\?[0-9]*" \
         | wc -l) -eq 0 ]] ; then
    echo "error in $0: ARGVALUE=$ARGVALUE is not a number" >&2
    exit 1
  fi
  case "$ARGUNIT" in
    H) NEWVALUE=$(echo "scale=8; $ARGVALUE * $HARTREE " | bc -l) 
       ;;
    RY) NEWVALUE=$(echo "scale=8; $ARGVALUE * $RYDBERG " | bc -l) 
       ;;
    EV) NEWVALUE=$(echo "scale=8; $ARGVALUE * $EV " | bc -l) 
       ;;
    KJ/MOL) NEWVALUE=$(echo "scale=8; $ARGVALUE * $KJBYMOL " | bc -l) 
       ;;
    KCAL/MOL) NEWVALUE=$(echo "scale=8; $ARGVALUE * $KCALBYMOL " | bc -l) 
       ;;
    ANGSTROM) NEWVALUE=$(echo "scale=8; $ARGVALUE * $ANGSTROM " | bc -l) 
       ;;
    CUBANGSTROM) NEWVALUE=$(echo "scale=8; $ARGVALUE * $CUBANGSTROM " | bc -l) 
       ;;
    ABOHR) NEWVALUE=$(echo "scale=8; $ARGVALUE * $ABOHR " | bc -l) 
       ;;
    CUBABOHR) NEWVALUE=$(echo "scale=8; $ARGVALUE * $CUBABOHR " | bc -l) 
       ;;
    *)
      echo "error in $0; (in function toatomicunits)" >&2
      echo "ARGUNIT=$ARGUNIT not recognized" >&2
      exit 1  
      ;;
  esac
  case "$NEWUNIT" in
    H) NEWVALUE=$(echo "scale=8; $NEWVALUE / $HARTREE " | bc -l) 
       ;;
    RY) NEWVALUE=$(echo "scale=8; $NEWVALUE / $RYDBERG " | bc -l) 
       ;;
    EV) NEWVALUE=$(echo "scale=8; $NEWVALUE / $EV " | bc -l) 
       ;;
    KJ/MOL) NEWVALUE=$(echo "scale=8; $NEWVALUE / $KJBYMOL " | bc -l) 
       ;;
    KCAL/MOL) NEWVALUE=$(echo "scale=8; $NEWVALUE / $KCALBYMOL " | bc -l) 
       ;;
    ANGSTROM) NEWVALUE=$(echo "scale=8; $NEWVALUE / $ANGSTROM " | bc -l) 
       ;;
    CUBANGSTROM) NEWVALUE=$(echo "scale=8; $NEWVALUE / $CUBANGSTROM " | bc -l) 
       ;;
    ABOHR) NEWVALUE=$(echo "scale=8; $NEWVALUE / $ABOHR " | bc -l) 
       ;;
    CUBABOHR) NEWVALUE=$(echo "scale=8; $NEWVALUE / $CUBABOHR " | bc -l) 
       ;;
    *)
      echo "error in $0; (in function toatomicunits)" >&2
      echo "NEWUNIT=$NEWUNIT not recognized" >&2
      exit 1  
      ;;
  esac
  #
  echo "$NEWVALUE $NEWUNIT"
}
################################################################################
#  resolve input arguments
#  $TN is "T" fort text output and "N" for number output
################################################################################
export TN=T
export WHATID="etot"
export VERBOSE=no
export UNIT=""
while getopts :hw:tnp:vu: OPT ; do
  case $OPT in 
    w) WHATID="$OPTARG" 
       case $WHATID in
         etot|gap|homo|lumo|efermi|volume|nkpt|kdiv|epw|epwrho) 
         ;; 
         *)
           echo "error in $0: illegal argument $WHATID of option -w" >&2
           exit 1
         ;;
       esac
       ;;
    n) TN=N             ;;
    t) TN=T             ;;
    u) UNIT=$OPTARG     
       UNIT=$(echo $UNIT | tr '[:upper:]' '[:lower:]') #convert to lower case
       case $UNIT in 
         h)        UNIT="H" ;;
         ry)       UNIT="RY" ;;
         ev)       UNIT="EV" ;;
         kj/mol)   UNIT="KJ/MOL" ;;
         kcal/mol) UNIT="KCAL/MOL" ;;
         angstrom^3|cubangstrom|aa^3) UNIT="CUBANGSTROM" ;;
         abohr^3|cubabohr|a0^3) UNIT="CUBABOHR" ;;
         *) echo "error in $0: argument $OPTARG for -u not recognized" >&2
             exit 1 ;;
       esac
       ;;
    v) VERBOSE=yes      ;;
    h) echo -e $USAGE ; exit 0  ;;
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
ROOT=$1
#
# report input variables -------------------------------------------------------
#
if [[ $VERBOSE = yes ]] ; then
  echo WHATID..................: $WHATID
  echo PROJECT NAME............: $ROOT
  echo UNIT....................: $UNIT
  echo VERBOSE.................: $VERBOSE
  echo text or numerical output: $TN
fi
#
################################################################################
#  check input
################################################################################
#
# is rootname defined ? --------------------------------------------------------
#
if [[ -z $ROOT ]] ; then
  echo error: missing argument. specify root name of the paw_project >&2
  exit 1
fi
#
#  is the output unit been specified? ------------------------------------------
#
if [[ -z $UNIT ]] ; then
  case $WHATID in 
    nkpt) ;;
    kdiv) ;;
    *) echo error: missing argument -u. specify UNIT  >&2 ; exit 1 ;;
  esac
fi
#
# Does protocol file exist?
#
if [[ ! -e $ROOT.prot ]] ; then
  echo "error in $0: protocol file does not exist" >&2
  echo "specified protocol file: $ROOT.prot" >&2
  exit 1
fi
#
# Is the unit consistent with quantity?
#
case $WHATID in
  etot|gap|homo|lumo|efermi|epw|epwrho)
    case $UNIT in
      H|EV|KJ/MOL|KCAL/MOL|RY) ;;
      *)
        echo "error in $0: unit inconsistency" >&2
        echo "energy variable $WHATID inconsistent with unit $UNIT" >&2
        exit 1
        ;;
    esac
    ;;
  volume) 
    case $UNIT in
      ABOHR^3|CUBABOHR|ANGSTROM^3|CUBANGSTROM) ;;
      *)
        echo "error in $0: unit inconsistency" >&2
        echo "volume variable $WHATID inconsistent with unit $UNIT" >&2
        exit 1
        ;;
    esac
    ;;
  nkpt|kdiv)
    if [[ -n $UNIT ]] ; then
      echo "error in $0: unit inconsistency" >&2
      echo "volume variable $WHATID inconsistent with unit $UNIT" >&2
      exit 1
    fi
    ;;
esac
#
################################################################################
#  collect data in a.u.
################################################################################
#
#--------------  total energy --------------------------------------------------
#
if [[ "$WHATID" = "etot" ]]  ; then 
  WHAT="TOTAL ENERGY"
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
#
#--------------  fermi level ---------------------------------------------------
#
elif [[ "$WHATID" = "efermi" ]]  ; then 
  WHAT="CHEMICAL POTENTIAL" 
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  if [[ -z  "$RESULT" ]] ; then # fall-back for insulators
    export WHAT1="HOMO-ENERGY" 
    export RESULT1=$(grep "$WHAT1" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
    RESULT1=${RESULT1%EV*}
    export WHAT2="LUMO-ENERGY" 
    export RESULT2=$(grep "$WHAT2" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
    RESULT2=${RESULT2%EV*}
    RESULT=$(echo " scale=5; ( $RESULT1 + $RESULT2 ) / 2 " | bc)
    RESULT="$RESULT EV"  
    if [[ -z $RESULT1 ]] ; then RESULT='';  fi
    if [[ -z $RESULT2 ]] ; then RESULT='';  fi
  fi
#
#--------------  band gap ------------------------------------------------------
#
elif [[ "$WHATID" = "gap" ]] ; then
  WHAT="ABSOLUTE GAP"
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  if [[ -z $RESULT ]] ; then  
    echo "error in $0: variable occupations: gap is unknown."  >&2
    exit 1
  fi
#
#--------------  homo  ---------------------------------------------------------
#
elif [[ "$WHATID" = "homo" ]] ; then
  export WHAT="HOMO-ENERGY"  
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  if [[ -z $RESULT ]] ; then  # fall back for metals
    export WHAT="CHEMICAL POTENTIAL"  
    export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  fi
#
#--------------  lumo  ---------------------------------------------------------
#
elif [[ "$WHATID" = "lumo" ]] ; then
  export WHAT="LUMO-ENERGY" 
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  if [[ -z $RESULT ]] ; then  # fall back for metals
    export WHAT="CHEMICAL POTENTIAL"  
    export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  fi
#
#--------------  volume  -------------------------------------------------------
#
elif [[ "$WHATID" = "volume" ]] ; then
  ##  cannot grep for "[". therefore enclose in "[ ]", i.e. [[]
  WHAT='[[]ANGSTROM[]]'
  export T11=$(grep "T1$WHAT" $ROOT.prot | tail -n1 | awk '{print $2}')
  export T12=$(grep "T1$WHAT" $ROOT.prot | tail -n1 | awk '{print $3}')
  export T13=$(grep "T1$WHAT" $ROOT.prot | tail -n1 | awk '{print $4}')
  export T21=$(grep "T2$WHAT" $ROOT.prot | tail -n1 | awk '{print $2}')
  export T22=$(grep "T2$WHAT" $ROOT.prot | tail -n1 | awk '{print $3}')
  export T23=$(grep "T2$WHAT" $ROOT.prot | tail -n1 | awk '{print $4}')
  export T31=$(grep "T3$WHAT" $ROOT.prot | tail -n1 | awk '{print $2}')
  export T32=$(grep "T3$WHAT" $ROOT.prot | tail -n1 | awk '{print $3}')
  export T33=$(grep "T3$WHAT" $ROOT.prot | tail -n1 | awk '{print $4}')

  # text wheter info exits
  if [[ -z $T11 ]] ; then
    echo "error in $0: no lattice constants in prot file" >&2
    echo "ROOT.prot= $ROOT.prot"
    exit 1
  fi


  DETSTRING=" $T11 * ( $T22 * $T33 - $T23 * $T32 ) \
            + $T12 * ( $T23 * $T31 - $T21 * $T33 ) \
            + $T13 * ( $T21 * $T32 - $T22 * $T31 ) "
  export RESULT=$(echo " scale=8;  ${DETSTRING} " | bc -l )
  RESULT="$RESULT CUBANGSTROM"
#
#--------------  nkpt  ---------------------------------------------------------
#
elif [[ "$WHATID" = "nkpt" ]] ; then
  export WHAT="NUMBER OF K-POINTS" 
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
#
#--------------  nkpt  ---------------------------------------------------------
#
elif [[ "$WHATID" = "kdiv" ]] ; then
  export WHAT="DIVISIONS FOR K-GRID" 
  export NK1=$(grep "$kdiv(1)" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  export NK2=$(grep "$kdiv(2)" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  export NK3=$(grep "$kdiv(3)" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
  export RESULT="${NK1} ${NK2} ${NK3}"
#
#--------------  EPW(psi)-------------------------------------------------------
#
elif [[ "$WHATID" = "epw" ]] ; then
  export WHAT="WAVEFUNCTION PLANE WAVE CUTOFF" 
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
#
#--------------  EPW(rho)-------------------------------------------------------
#
elif [[ "$WHATID" = "epwrho" ]] ; then
  export WHAT="DENSITY PLANE WAVE CUTOFF" 
  export RESULT=$(grep "$WHAT" $ROOT.prot | tail -n1 | awk -F: '{print $2}')
#
else
  echo "error in $0: $WHAT is not a permitted value for WHAT"  >&2
  exit 1
fi
#  test if result is empty
if [[ -z  "$RESULT" ]] ; then
  echo "error in $0: no result found in $ROOT.prot" >&2
  exit 1
fi
#
################################################################################
#  convert units
################################################################################
case $WHATID in
 nkpt|kdiv) ;;                            # no unit 
 *) RESULT=$(convertunit $UNIT "$RESULT") # internal function convertunit 
   ;;
esac
#
################################################################################
#  switch between number and text output and edit result
################################################################################
if [[ $TN = N ]] ; then 
  case $WHATID in
    kdiv) ;;
    *) RESULT=$(echo "$RESULT" | awk '{print $1}') ;;
  esac
else
  RESULT=$(echo "$WHATID = $RESULT for $ROOT")
fi
#
################################################################################
#  report result
################################################################################
echo "$RESULT"
exit 0
