#!/bin/bash
# this script changes the filenames of module file "FILE.mod" into all lowercase.
#
# Usage:  
#      lowercasemodnames.sh <filein >fileout
#
IFS=""
## Die Datei aus $1 zeilenweise auslesen und in "line" hinterlegen
while read -r linein ;  do 
  # ${string#substring} Strips shortest match of $substring from front of $string.
  # ${string##substring} Strips longest match of $substring from front of $string.
  # ${string%substring} Strips shortest match of $substring from back of $string.
  # ${string%%substring} Strips longest match of $substring from back of $string.
  # ${#string} length of $string
  #
  # test if a module file is on this line
  X=`echo ${linein} | grep "\.mod"`
  if test ${#X} -eq 0  ; then
    # no module file; just copy input into output
    echo $linein
  else
    # referral to a module file 
    # assumes that line has the form of a dependency, i.e. FIELD1 : FIELD2
    FIELD1=${linein%:*}
    STRIP=${FIELD1%.mod*}
    if test ${#STRIP} -ne ${#FIELD1} ; then
      FIELD1=`echo $STRIP | tr ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz` 
      FIELD1=${FIELD1}.mod
    fi    
    FIELD2=${linein#*:}
    STRIP=${FIELD2%.mod*}
    if test ${#STRIP} -ne ${#FIELD2} ; then
      FIELD2=`echo $STRIP | tr ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz` 
      FIELD2=${FIELD2}.mod
    fi    
    lineout="${FIELD1} : ${FIELD2}"
    echo "$lineout"
  fi
done 
