#!/bin/bash
######################################################################
#
#  name:    paw_copy.sh
#
#  purpose: create a copy of a cp-paw project with a new root name
#
#  usage:  paw_copy(.sh) name1 name2
#    name1: root name of the existing cp-paw project
#    name2: root name of the new cp-paw project
# 
# copies all files beginning with string $1 into files beginning with
# $2 and with teh same ending
#
#######################################################################
for FILE in $1* ;do 
  cp -p "$FILE" "$2${FILE#$1}"
  RC=$?
  if [[ $RC -ne 0 ]] ; then 
     echo "error in $0" >&2
     echo "copying $FILE to $2${FILE#$1} failed" >&2
     exit 1
  fi
done
exit 0
