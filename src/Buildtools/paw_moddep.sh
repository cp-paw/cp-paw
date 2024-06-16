#!/bin/bash
################################################################################
##                                                                            ##
##    paw_moddep.sh                                                           ##
##                                                                            ##
##    Purpose:                                                                ##
##      outputs the make-dependencies of module files from fortran code       ##
##                                                                            ##
##    Usage:                                                                  ##
##      paw_moddep.sh options -f file  > dependencies                         ##
##                                                                            ##
##      file is a fortran code (Files not ending on .[Ff]* will be ignored)   ##
##                                                                            ##
##      the dependencies are written to standard output and can be passed     ##
##      into a make file                                                      ##
##                                                                            ##
##      dependencies is a file collecting the dependencies to be introduces   ##
##      into a make file                                                      ##
##                                                                            ##
##    Author: Peter Bloechl, Goslar 2024                                      ##
##                                                                            ##
################################################################################

#-------------------------------------------------------------------------------
#  set up help message $USAGE
#-------------------------------------------------------------------------------
export USAGE="\n"
USAGE="$USAGE Usage of $0:\n"
USAGE="$USAGE \n"
USAGE="$USAGE \t paw_moddep.sh options -f filename >dependencies\n"
USAGE="$USAGE \n"
USAGE="$USAGE extracts the module-related dependencies from a 
              fortran source file.\n"
USAGE="$USAGE The dependency list is written to standard out. \n"
USAGE="$USAGE The dependency list can be integrated into a make file. \n"
USAGE="$USAGE \n"
USAGE="$USAGE Files without a fortran extension e.g. f,f90,f08 etc are ignored.\n"
USAGE="$USAGE Dependencies on module files within the same file are omitted.\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options \n"
USAGE="$USAGE \t -f filename\t name of the fortran file to be parsed.\n"
USAGE="$USAGE \t -d directory\t directory for module and object files.\n"
USAGE="$USAGE \t -u \t\t uppercase module file names e.g. XXX.mod.\n"
USAGE="$USAGE \t\t\t The default are lowercase file names e.g. xxx.mod\n"
USAGE="$USAGE \t -v \t\t verbose\n"
USAGE="$USAGE \t -h \t\t prints this help message\n"

#-------------------------------------------------------------------------------
#  resolve options
#-------------------------------------------------------------------------------
export TUPPERMOD=""  # set to yes for XXX.mod rather than xxx.mod
export VERBOSE=""
export IN=""
while getopts :hvuf:d: OPT ; do
  case $OPT in
    f) IN=$OPTARG ;;
    d) OBJDIR=$OPTARG ;;
    u) TUPPERMOD="Y" ;;
    v) VERBOSE="Y" ;;
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
if [[ -z $IN ]] ; then 
   echo "error in $0: file not specified"                                    >&2
   exit 1
fi 
if [[ -z $OBJDIR ]] ; then 
   echo "error in $0: directory for module and object files not specified"   >&2
   exit 1
fi 
if [[ -n $VERBOSE ]] ; then
   echo INPUT FILE NAME = $IN
   echo directory for module and object files = $OBJDIR
fi

#-------------------------------------------------------------------------------
#-- skip nonfortran files                                                     --
#-------------------------------------------------------------------------------
X=${IN#*.}     # keep only the extension 
if [[ ! ${X} = [Ff][0-9][0-9] ]] ; then  
  if [[ -n $VERBOSE ]] ; then
     echo "INPUT FILE is not a fortran file"
     echo "extension $X"
  fi
  exit 0 
fi

export OBJ
#OBJ=${OBJDIR}/${IN##*/}
OBJ=${IN##*/}
OBJ=${OBJ%.f*}
OBJ=${OBJ%_d}.o    # allow also _d.f90 files

export TARGETS="${OBJ}"
export PREREQUISITES="${IN}"
#-------------------------------------------------------------------------------
#-- collect dependencies of module files on fortran files                     --
#-------------------------------------------------------------------------------
export TMPFILE=$(mktemp)
grep -i module $IN > ${TMPFILE}
LIST=""
while read -r FIRST SECOND REST ; do
  FIRST=$(echo "$FIRST" | tr '[:upper:]' '[:lower:]' )
  SECOND=$(echo "$SECOND" | tr '[:upper:]' '[:lower:]' )
  if [[ -n $VERBOSE ]] ; then
    echo FIRST SECOND REST=$FIRST $SECOND $REST
  fi
  if [[ "${FIRST}" = "module" ]] ; then
    SECOND=${SECOND%%!*}   # remove trailing "!"
     
    SKIP=""
    if [[ $SECOND == procedure ]] ; then 
       SKIP=true 
    fi

    if [[ -z $SKIP ]] ; then
      LIST="$LIST $SECOND "  # add to list of modules
      if [[ -n $TUPPERMOD ]] ; then 
        SECOND=$(echo "$SECOND" | tr '[:lower:]' '[:upper:]')
      fi
#      echo ${OBJDIR}/${SECOND}.mod : ${OBJ} >> /dev/stdout
#      TARGETS="${TARGETS} ${OBJDIR}/$SECOND.mod"
      TARGETS="${TARGETS} $SECOND.mod"
    fi
  fi
done <  "${TMPFILE}"

#-------------------------------------------------------------------------------
#-- collect dependencies of fortran files on module files                     --
#-------------------------------------------------------------------------------
grep -i use $IN > ${TMPFILE}
while read -r FIRST SECOND REST ; do
  FIRST=$(echo "$FIRST" | tr '[:upper:]' '[:lower:]' )
  SECOND=$(echo "$SECOND" | tr '[:upper:]' '[:lower:]' )
  # echo FIRST=$FIRST
  # echo SECOND=$SECOND
  # echo REST=$REST
  if [[ "${FIRST}" = "use" ]] ; then
    SECOND=${SECOND%%,*}   # remove trailing ", only:" 
    SECOND=${SECOND%%!*}   # remove trailing "!"

    SKIP=""
    # exclude dependencies on modules in the same file
    for X in ${LIST} ; do
       if [[ $SECOND == $X ]] ; then SKIP=true ; break ; fi
    done

    # exclude module files from libraries included via FCFLAGS etc.
    for X in xc_f03_lib_m.mod mpi_f08.mod ; do
       if [[ $SECOND == ${X%.mod} ]] ; then SKIP=true ; break ; fi
    done

    if [[ -z $SKIP ]] ; then
      LIST="$LIST $SECOND"
      if [[ -n $TUPPERMOD ]] ; then 
        SECOND=$(echo "$SECOND" | tr '[:lower:]' '[:upper:]')
      fi
#      echo  ${OBJ} : ${OBJDIR}/${SECOND}.mod >> /dev/stdout
      PREREQUISITES="${PREREQUISITES} ${SECOND}.mod "
    fi
  fi
done <  "${TMPFILE}"

#-------------------------------------------------------------------------------
#-- clean up                                                                  --
#-------------------------------------------------------------------------------
echo "${TARGETS} &: ${PREREQUISITES}" >> /dev/stdout
echo -e '\t $(FC) -c $(FCFLAGS)' -o${OBJ} ${IN} >> /dev/stdout
# the date of the module files is not updated if it needs no rebuild so that
# a the rule may be applied again.
echo -e '\t touch' ${TARGETS} >> /dev/stdout
echo -e "##" >> /dev/stdout

#-------------------------------------------------------------------------------
#-- clean up                                                                  --
#-------------------------------------------------------------------------------
rm ${TMPFILE}

