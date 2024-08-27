#!/bin/bash
###################################################################
##                                                               ##
##  name: paw_show                                               ##
##                                                               ##
##  purpose: pulls out columns from the paw protocol             ##
##                                                               ##
##  usage:                                                       ##
##    pawshow option --- rootname                                ##
##                                                               ##
##  options:                                                     ##
##     -h       help                                             ##
##     -e       static total energy                              ##
##     -t       atomic temperature                               ##
##     -c       constant energy (constant of motion)             ##
##     -f       fictitios kinetic energy of the wave functions   ##
##     -s       time steps rather than time in ps                ##
##     -u unit  select unit (h, ry, ev, kj/mol, kcal/mol)        ##
##     -a type: friction                                         ##
##         type=r   friction acting on nuclei                    ##
##         type=p   friction acting on wave functions            ##
##     -o file: produce output file rather than opening xmgrace  ##
##              extension (eps,png) specifies graphics format    ##
##              other (or no) extensions produce data for grace  ##
##              execute with xmgrace -nxy name                   ##
##              if option -p is absent, xmgrace is started       ##
##     -b batch:add settings from batch file to xmgrace settings ##
##                                                               ##
##  dpendency:                                                   ##
##    uses xmgr as plotting routine                              ##
##                                                               ##
###################################################################
USAGE="
Usage: paw_show options --- rootname\n
1)rootname  is paw protocol file name without the .prot ending 
2)option can be one or more of the following 
\t -h: give info about use and exit 
\t -s: use time steps rather than time in picoseconds
\t -u unit: select energy unit (h, ry, ev, kj/mol, kcal/mol)      
\t -e: plot static total energy E_tot in Hartree (or unit)
\t -c: plot conserved energy E_c in Hartree (or unit)
\t -f: plot fictitious kinetic wave function energy in Hartree (or unit)
\t -t: plot temperature T in Kelvin
\t -a type: plot friction 
\t\t     type=p: for wave functions a_p
\t\t     type=r: for atoms atoms a_R 
\t -o filename: print to file rather than starting xmgrace. 
\t\t     extensions eps,png produce graphics files
\t\t     others or no extensions produce xmgrace input file
\t\t     process with xmgrace -nxy filename
\t -b batchfile: attach batch file to xmgrace settings
\n
Purpose: plot time-dependent data from protocol file
\n"
# use the following for qm-mm calculations:\n
# \t -qe: plot static total energy of the environmet\n
# \t -qc: plot conserved energy of the environment\n
# \t -qt: plot temperature in kelvin of the environment\n
# \t -tqt: both temperatures\n
# \t -cqc: both conserved energies\n
#===================================================================
#  resolve argument list                                          ==
#===================================================================
export XLABEL='t[ps]' # default
export UNIT='H'       # default
export EUNIT=1.       # default
export YLABEL=""     
declare -i ICOLOR=1
OPTIND=1
AWKSTRING=""
export PRINTFILE=""
export BATCHFILE=""
#
# xmgrace typesetting:
#   \R{n} switch to color n
#   \s subscript
#   \N normal font
#   \f{SYMBOL} switch to greek font
#   \f{} switch to normal font
#
while getopts :hetfca:su:o:b: OPT ; do
  case $OPT in 
    t) AWKSTRING="$AWKSTRING \" \" \$4 " # temperature (Kelvin)
       YLABEL="${YLABEL}\R{$ICOLOR}T[K]; "
       ICOLOR=$ICOLOR+1
       ;;
    f) AWKSTRING="$AWKSTRING \" \" \$5/$EUNIT "  # fict. kinetic energy
       YLABEL="${YLABEL}\R{$ICOLOR}E\s\f{Symbol}y\f{}\N[$UNIT]; "
       ICOLOR=$ICOLOR+1
       ;;
    e) AWKSTRING="$AWKSTRING \" \" \$6/$EUNIT " # static total energy
       YLABEL="${YLABEL}\R{$ICOLOR}E\stot\N[$UNIT]; "
       ICOLOR=$ICOLOR+1
       ;;
    c) AWKSTRING="$AWKSTRING \" \" \$7/$EUNIT " # conserved energy
       YLABEL="${YLABEL}\R{$ICOLOR}E\sc\N[$UNIT]; "
       ICOLOR=$ICOLOR+1
       ;;
    a) #friction
       case $OPTARG in
         p) AWKSTRING="$AWKSTRING \" \" \$8"  # wave functions
            YLABEL="${YLABEL}\R{$ICOLOR}a\s\f{Symbol}y\f{}\N; "
            ICOLOR=$ICOLOR+1
            ;;
         r) AWKSTRING="$AWKSTRING \" \" \$9"  # atomic positions
            YLABEL="${YLABEL}\R{$ICOLOR}a\sR\N; "
            ICOLOR=$ICOLOR+1
           ;;
         *) echo "error in $0: argument $OPTARG not recognized" >&2
            exit 1
            ;;
       esac
       ;;
    s) XLABEL="time steps" ;; # time coordinate 
    u) UNIT=$OPTARG 
       UNIT=$(echo $UNIT | tr '[:upper:]' '[:lower:]') #convert to lower case
       case $UNIT in
         h)        EUNIT=1.; UNIT="H" ;;
         ry)       EUNIT=0.5; UNIT="Ry" ;;
         ev)       EUNIT=$(echo "scale=8; 1. / 27.211 " | bc -l) 
                   UNIT="eV"
                   ;;
         kj/mol)   EUNIT=$(echo "scale=8; 1. / 2625.500223430069 " | bc -l) 
                   UNIT="kJ/mol"
                   ;;
         kcal/mol)EUNIT=$(echo "scale=8; 1. / 627.5096463920391 " | bc -l) 
                   UNIT="kcal/mol"
                   ;;
         *) echo "error in $0: argument $OPTARG for -u not recognized" >&2
             exit 1 ;;
       esac
      ;;
    o) # specify output file
       PRINTFILE=$OPTARG
       ;;
    b) # specify batchfile to add to parameter file for xmgrace
       BATCHFILE=$OPTARG
       ;;
    h) 
       echo -e "$USAGE"
       exit 0
       ;;
    \?) # unknown option
       echo "error in $0: invalid option -$OPTARG" >&2
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
#
ROOT=$1
if [[ -z $ROOT ]] ; then
  echo "error in $0: argument ROOT not specified" >&2
  exit 1
fi

# second argument (output file) is outdated. Fallback
if [[ $# > 1 ]] ; then
   echo "Warning: outdated option" >&2
   echo "second argument is treated as output file" >&2
   if [[ -n ${PRINTFILE} ]] ; then
     echo "error in $0" >&2
     echo "printfile set and second argument present" >&2
     exit 1
   fi
   export PRINTFILE=$2
fi
#
#===================================================================
#  make a temporary file
#===================================================================
# check if environment variable TMPDIR is set
if [[ -z ${TMPDIR} ]] ; then TMPDIR=/tmp ; fi
if [[ ! -d ${TMPDIR} ]] ; then
  echo "error in $0: temp directory does not exist"
  echo "TMPDIR=$TMPDIR"
  exit 1
fi
TMPFILE=$(mktemp ${TMPDIR}/pawshowfile.XXXXXX)
RC=$?
if [[ $RC -ne 0 ]]; then
  echo "error in $0: Can't create temp file, exiting..."
  exit 1
fi
PARMFILE=$(mktemp ${TMPDIR}/pawshowfile.XXXXXX)
RC=$?
if [[ $RC -ne 0 ]]; then
  echo "error in $0: Can't create temp file, exiting..."
  exit 1
fi
#
#===================================================================
#  select time in ps or time steps
#===================================================================
if [[ "$XLABEL" = "t[ps]" ]] ; then
  AWKSTRING="{print  \$3 $AWKSTRING}"
else
  AWKSTRING="{print  \$2 $AWKSTRING}"
fi
#
#===================================================================
#  extract data from protocol file                                ==
#===================================================================
grep "\!>" ${ROOT}.prot | sed 's/!>/!> /1' \
                       | sed 's/-/ -/g' \
                       | gawk -v CONVFMT=%.17g "$AWKSTRING" \
                       > ${TMPFILE}
#
#===================================================================
#  construct parameter file
#===================================================================
>$PARMFILE
echo 'map font 0 to "Helvetica", "Helvetica" '       >> $PARMFILE
echo 'default linewidth 3.0'                         >> $PARMFILE
echo 'page background fill off'                      >> $PARMFILE
echo 'frame background pattern 1'                    >> $PARMFILE
echo 'frame background color 0'                      >> $PARMFILE
echo 'default sformat "%.17g"'                        >> $PARMFILE
echo "with g0"                                       >> $PARMFILE
echo "  view 0.2,0.15,1.2,0.95"                      >> $PARMFILE
echo "  legend on"                                   >> $PARMFILE
if [[ "$XLABEL" = "t[ps]" ]] ; then
  echo '  xaxis LABEL "t[PS]"'                       >> $PARMFILE
else
  echo '  xaxis LABEL "time steps" '                 >> $PARMFILE
fi
echo "  yaxis ticklabel prec 14"                     >> $PARMFILE
YLABEL=${YLABEL%;*}   # remove trailin "; "
echo "  yaxis LABEL \"${YLABEL}\""                   >> $PARMFILE
#===================================================================
#  append instructions from external batch file to the parameterfile
#  this should overwrite previous settings and allows the user to 
#  automate the plot appearance
#===================================================================
if [[ -n $BATCHFILE ]] ; then
  cat $BATCHFILE >> $PARMFILE
fi
#
#===================================================================
#  print or present via xmgrace
#===================================================================
if [[ -n ${PRINTFILE} ]] ; then
    # gracebat -nosafe -hdevice EPS -${GRACEINCLUDE} $PARMFILE \
    #          -hardcopy -nxy ${TMPFILE} -printfile ${PRINTFILE}.eps
  case ${PRINTFILE##*.} in
    eps)
      export DEVICE="EPS"
      ;;
    png)
      export DEVICE="PNG"
      ;;
    *)  # write xmgrace input data to file
        # execute with "xmgrace -nxy file"
      export DEVICE=""
      cp ${TMPFILE} $PRINTFILE
      ;;
  esac
  if [[ -n $DEVICE ]] ; then
    gracebat -nosafe -hdevice ${DEVICE} -batch $PARMFILE \
             -hardcopy -nxy ${TMPFILE} -printfile ${PRINTFILE}
  fi
else     #-- open xmgrace if no printfile specified
  xmgrace -free -noask -batch $PARMFILE -nxy ${TMPFILE}
fi
rm ${TMPFILE}
rm ${PARMFILE}
exit 0




