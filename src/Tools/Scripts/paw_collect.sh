#!/bin/bash
#################################################################
#
#  name: paw_collect.sh 
#
#  Usage: paw_collect(.sh)
#
#  purpose: collects total energies from CP-PAW projects
#
# searches through all *.prot files in the current directory, 
# as well as up to four subdirectories down
# and lists the name of the project and the last total energy
# written in the extended printout. 
#################################################################
export USAGE="Usage:
\t     paw_collect options 
Options:
\t  -h: print this help message and exit
Prints the total energies of all paw_projects in the current directory 
and five directories down. The unit is Hartree.
\n"

while getopts :h OPT ; do
  case $OPT in
    h)
      echo -e "$USAGE"
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

for FILE in         *.prot \
                  */*.prot \
                */*/*.prot \
              */*/*/*.prot \
            */*/*/*/*.prot ; do
  if [[ -f $FILE ]] ; then
    ROOT=${FILE%.prot}
    ENERGY=$(paw_get.sh -nw etot -u H -- $ROOT)
    echo -e "$ROOT $ENERGY"
  fi
done
