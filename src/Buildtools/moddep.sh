#!/bin/bash
################################################################################
##                                                                            ##
##    moddep.sh                                                               ##
##                                                                            ##
##    Purpose:                                                                ##
##      outputs the make-dependencies of module files from fortran code       ##
##                                                                            ##
##    Usage:                                                                  ##
##      moddep.sh file                                                        ##
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
export IN=$1
#-------------------------------------------------------------------------------
#-- skip nonfortran files                                                     --
#-------------------------------------------------------------------------------
X=${IN#*.}     # keep only the extension 
if [[ ! ${X} = [Ff][0-9][0-9] ]] ; then  exit 0 ; fi

#-------------------------------------------------------------------------------
#-- collect dependencies of module files on fortran files                     --
#-------------------------------------------------------------------------------
export TMPFILE=$(mktemp)
grep -i module $IN > ${TMPFILE}
LIST=""
while read -r FIRST SECOND REST ; do
  FIRST=$(echo "$FIRST" | tr '[:upper:]' '[:lower:]' )
  # echo FIRST=$FIRST
  # echo SECOND=$SECOND
  # echo REST=$REST
  if [[ "${FIRST}" = "module" ]] ; then
    SECOND=${SECOND%%!*}   # remove trailing "!"
     
    SKIP=""
    if [[ $SECOND == [Pp][Rr][Oo][Cc][Ee][Dd][Uu][Rr][Ee] ]] ; then 
       SKIP=true 
    fi

    if [[ -z $SKIP ]] ; then
      LIST="$LIST $SECOND "  # add to list of modules
      echo $SECOND.mod : $IN >> /dev/stdout
    fi
  fi
done <  "${TMPFILE}"

#-------------------------------------------------------------------------------
#-- collect dependencies of fortran files on module files                     --
#-------------------------------------------------------------------------------
grep -i use $IN > ${TMPFILE}
while read -r FIRST SECOND REST ; do
  FIRST=$(echo "$FIRST" | tr '[:upper:]' '[:lower:]' )
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
    if [[ -z $SKIP ]] ; then
      LIST="$LIST $SECOND"
      echo  $IN : $SECOND.mod >> /dev/stdout
    fi
  fi
done <  "${TMPFILE}"

#-------------------------------------------------------------------------------
#-- clean up                                                                  --
#-------------------------------------------------------------------------------
rm ${TMPFILE}

