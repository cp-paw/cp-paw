#!/bin/bash
################################################################################
#
#         FILE: paw_resolve.sh
#
#        USAGE:  paw_resolve(.sh) options infile > outfile
#
#      OPTIONS: -h, -p prefix, -r rule, -f frule
#
#  DESCRIPTION: 
#
#  parses the input file for occurances of the type @string@ and
#  replaces them in one of two possible ways:
#
#  1) if string is a keyword that occurs in an explicit rule supplied
#  to the argument with the -r option, @string@ is replaced by the
#  specified replacement.
#
#  2) if a string is a keyword that occurs in an frule (file-rule)
#  supplied to the argument with the -f option, @string@ is replaced
#  by the contents of the specified file. CAUTION: the complete line
#  containing $string$ is removed in the process
#
#  3) if string is not specified as a rule, it is considered the name
#  of a file, after adding the prefix supplied with the option -p.
#  The line containing @string@ is removed and replaced by the
#  contents of the file ${PRE}string.  The prefix can be used to
#  specify a directory that holds the insertions.
#
#  Remark: @string@ must not contain any spaces. Two neighboring
#  occurances must be separated by a space or a newline.
#
#  Remark: colons are used instead of slashes as separators in sed to allow
#  replacements with path names.
#
#  Remark: string replacement rules (-r) precede file replacement rules (-f).
#  i.e. they do not act on inserted files
#
#  Example:
#
#    paw_resolve.sh -r id1=rep1 -rid2=rep2 -pInsertions/ infile >outfile
#
#       AUTHOR: Peter Bloechl; peter.bloechl@tu-clausthal.de
#      CREATED: 14. Dec. 2013
#
################################################################################
#
# description of usage
#
USAGE="Usage of paw_resolve.sh:\n"
USAGE="$USAGE \n"
USAGE="$USAGE \t paw_resolve(.sh) options\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -i input file (mandatory)\n"
USAGE="$USAGE \t -o output file (default: stdout)\n"
USAGE="$USAGE \t -p prefix: prefix to be added to the filename in @file@\n"
USAGE="$USAGE \t -r id=value: rule to convert @id@ into value\n"
USAGE="$USAGE \t -f id=file: rule to replace lines containing  @id@ "
USAGE="$USAGE       by the contents of file \n"
USAGE="$USAGE \t -h this help message\n"
USAGE="$USAGE \t -v verbose\n"
USAGE="$USAGE \n"
USAGE="$USAGE \t Caution: specify output file with -o to when using -v or -h,"
USAGE="$USAGE \n\t\t to avoid piping undesired information into the outfile"
#-------------------------------------------------------------------------------
# resolve argument list
#-------------------------------------------------------------------------------
INFILE=""      # input file
OUTFILE=""     # output file
PRE=""         # prefix for file insertions (such as a path)
declare -a RULES=( ) # string replacement rules  A=B: @A@ -> B as array
    # declare -a RULES=(a b c d) create array RULES with four elements a b c d
    # ${RULES[3]}   Third array element 
    # ${RULES[3]}   Third array element 
    # ${RULES[*]}   all array elements
    # ${#RULES[*]}  number of array elements  
    # ${!RULES[*]}  array of all array indices
NRULE=0
FRULES="" # file insertion rules
INFILE=""
VERBOSE=no
while getopts :hi:o:p:r:f:v OPT ; do
  case $OPT  in
    i)
      INFILE="$OPTARG"
      ;;
    o)
      OUTFILE="$OPTARG"
      ;;
    p)
      PRE="$OPTARG"
      ;;
    r)
      let I=${#RULES[*]}+1 # number of array elements +1
      RULES[$I]=$OPTARG    # add $OPTARG as next element to array RULES
      ;;
    f)
      FRULES="$FRULES $OPTARG"
      ;;
    v)
      VERBOSE=yes
      ;;
    h)
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
#
# check if mandatory arguments are provided
#
if [[ -z "${INFILE}" ]] ; then
  echo "error in $0: mandatory argument -i infile missing" >&2
  exit 1
fi
#
#   report setting of parameters
#
if [[ $VERBOSE = yes ]] ; then
 echo input file : $INFILE
 if [[ -n "$PRE" ]]    ; then echo "prefix: $PRE"        ; fi
 for I in ${!RULES[*]} ; do
    echo RULE $I:${RULES[$I]}
 done
 if [[ -n "$FRULES" ]] ; then echo "file rules: $FRULES" ; fi
fi


if [ ! -e ${INFILE} ] ; then 
  echo "error in $0: file of file rule $X does not exist" >&2
     exit 1
   fi

#
#-------------------------------------------------------------------------------
# explicit string replacement rules
#-------------------------------------------------------------------------------
CHS0=""
CHS1=""
CHS2=""
for I in ${!RULES[*]} ; do    # step through array indices
   A=${RULES[$I]%=*}
   A=$(echo $A | sed 's/^ *//') #remove preceeding spaces
   A=$(echo $A | sed 's/ *$//') #remove trailing spaces
   B=${RULES[$I]##*=}
   CHS0="$CHS0-e s:@$A@:\"$B\":g "
done
if [[ $VERBOSE = yes ]] ; then
  echo -e "explicit string replacement rules: CHS0=$CHS0"
fi
#-------------------------------------------------------------------------------
# explicit file replacement rules
#-------------------------------------------------------------------------------
for X in $FRULES ; do
   A=${X%=*}
   B=${X##*=}
   if [ ! -e ${PRE}$B ] ; then 
     echo "error in $0: file of file rule $X does not exist" >&2
     exit 1
   fi
   CHS0="$CHS0 -e/@$A@/r${PRE}${B} "
   CHS2="$CHS2 -e/@$A@/d " 
done
#
#-------------------------------------------------------------------------------
# collect all lines containing @*@ in infile ($INFILE) and obtain the file
# names in @file@ as $LIST
#-------------------------------------------------------------------------------
VAR=$(sed -n -e '/@*@/p' < $INFILE)
# collect the filennames enclosed in @file@ in LIST
LIST=""
for X in $VAR; do
  Y=$X
  X=${X%@*}    # strip trailing '@' and following stuff
  X=${X#*@}    # strip preceeding '@' and everything before
  if [ $X = $Y ] ; then X="" ; fi  # discard expressions without @
  #
  #   check if $X corresponds to a direct replacement (-r). otherwise clear it
  for I in ${!RULES[*]} ; do
    A=${RULES[$I]%=*}
    A=$(echo $A | sed 's/^ *//') #remove preceeding spaces
    A=$(echo $A | sed 's/ *$//') #remove trailing spaces
    if [ "$A" = "$X" ] ; then X="" ; fi
  done
  for Y in $FRULES ; do
    A=${Y%=*}
    if [ "$A" = "$X" ] ; then X="" ; fi
  done
  #
  # avoid double counting
  for Y in $LIST; do
    if [ "$Y" = "$X" ] ; then X="" ; fi
  done
  #
  #
  #
  if [[ ! -z "$X" ]] ; then
    # check if there $X is the name of an existing file
    if [[ ! -e ${PRE}${X} ]] ; then
      echo "error in $0: no rule or file to replace @$X@" >&2
      exit 1
    fi
    # check if the lines containing file replacements are alone on one line
    LINES=$(sed -n -e ":@$X@:=" $INFILE)
    for Y in $LINES ; do
      LINE=$(sed -n -e ${Y}p < $INFILE)
      for Z in $LINE ; do
        if [[ "@$X@" != "$Z" ]] ; then
          echo "error in $0: file replacement @$X@ is not alone in the line" >&2
          echo "LINENO=$Y LINE=$LINE" >&2
          exit 1
        fi
      done
    done
    LIST="$LIST $X"
  fi
done
#-------------------------------------------------------------------------------
# replace @file@ with file and then remove the corresponding entries
# @file@ removal follows all replacements to deal properly with
# multiple occurances.  The $IN is processed and the result is written
# to standard out.
#-------------------------------------------------------------------------------
for X in $LIST; do
  CHS1="${CHS1}-e:@$X@:r${X}"
done
#
for X in $LIST; do
   CHS2="$CHS2 -e :@$X@:d " 
done
#-------------------------------------------------------------------------------
# compose and execute command
#-------------------------------------------------------------------------------
 # echo CHS0=$CHS0 ${CHS0[0]}
 # echo CHS1=$CHS1
 # echo CHS2=$CHS2
CMD="cat ${INFILE}"
if [[ -n "${CHS0}" ]] ; then CMD="$CMD | sed ${CHS0}" ; fi # file replacement
if [[ -n "${CHS1}" ]] ; then CMD="$CMD | sed ${CHS1}" ; fi
if [[ -n "${CHS2}" ]] ; then CMD="$CMD | sed ${CHS2}" ; fi
if [[ ! -z "${OUTFILE}" ]] ; then
  CMD="$CMD > ${OUTFILE}"
fi
if [[ $VERBOSE = yes ]] ; then echo $CMD; fi
eval ${CMD}
#-------------------------------------------------------------------------------
exit 0


