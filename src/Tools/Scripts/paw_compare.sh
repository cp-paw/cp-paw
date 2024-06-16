#!/bin/bash
################################################################################
#
#  paw_compare compares the fortran (.f90) files in two specified directories.
#  
#  1) It lists files that are present in only one of the two directories and
#  2) it identifies files that differ, but exist is both directories.
#  3) for the latter it will print the subroutine names that are affected
#  
#  for further information, I recommend the open source tool "opendiff"
#  
################################################################################
export THISDIR=$(pwd)
if [[ $# != 2 ]] ; then
  echo "error in $0: must have two arguments"
  echo "number of arguments is $#"
  exit 1
fi

export DIR1=$1
export DIR2=$2
#-- convert relative path into absolute path
cd ${THISDIR}
cd $DIR1; DIR1=$(pwd)
cd ${THISDIR}
cd $DIR2; DIR2=$(pwd)
#-- create temporary work directory
export WORKDIR=$(mktemp -d)

export WORKDIR="${THISDIR}/Here"
mkdir $WORKDIR
#
# list all .f90 files in .dir1 and .dir2
#
cd $DIR1; ls -1 *.f90 > $WORKDIR/.dir1
cd $DIR2; ls -1 *.f90 > $WORKDIR/.dir2
#
# check for differences in .dir1 and .dir2 and report
# inconsistencies concerning the filenames
#
cd $WORKDIR
diff -bBi --minimal -s  .dir1 .dir2 > .tmp
echo "----------- FILE REPORT ----------------"
echo
echo "==> The following *.f90 files of directory $DIR1 are not in $DIR2:"
echo
grep '<' .tmp | cut -d'<' -f2
echo
echo "==> The following *.f90 files of directory $DIR2 are not in $DIR1:"
echo
grep '>' .tmp | cut -d'>' -f2
echo

#
# create a list ($LIST) which just holds the common filenames
# of both directories!
#
LIST0=$(grep '<' .tmp | cut -d'<' -f2)
if [[ "$LIST0" == " " ]] ; then
  for i in $LIST0; do
    var='-e'"$i"
    LIST1="$LIST1 $var"
  done
  LIST=$(grep -v $LIST1 .dir1)
else
  LIST=$(cat .dir1)
fi

#
# LOOP over all common files.
#
for i in $LIST; do
  # The comment (!) lines are not evaluated
  #cat $DIR1/$i | grep -v ^! > .file1
  #cat $DIR2/$i | grep -v ^! > .file2
  cat $DIR1/$i  > .file1
  cat $DIR2/$i  > .file2
  #
  # check for differences 
  #
  diff -bBi --minimal -s .file1  .file2 > .tmp
  # 
  # if the output has more than one line (-> differences),
  # these the corresponding subroutine is searched
  #
  if [[ $(wc .tmp | awk '{print $1}') != "1" ]] ; then
    echo 
    echo "==> DIFFERENCES IN $i" 
    echo
    DIFFLIST=$(grep ^[0-9] .tmp  | cut -d'd' -f1 | cut -d'c' -f1  | cut -d',' -f1 | cut -d'a' -f1)
    SUBROUTINELINES=$(grep -n -eSUBROUTINE -eMODULE.*MODULE .file1 | grep -v 'END'| grep -v USE | cut -d':' -f1)

    k=1
    for j in $DIFFLIST; do
      FOUND="F"
      for k in $SUBROUTINELINES; do
        if test "$k" -gt "$j" -a "$FOUND" = "F"; then
          NAME=$(head -$kold .file1 | tail -1 | cut -d'(' -f1)
          if [[ "$NAME" != "$NAMEOLD" ]] ; then
            echo $NAME 
#	    echo line $j 
            NAMEOLD=$NAME
	  fi
	  echo -e "line $j: $(sed -n -e ${j}p ${DIR1}/$i)"
          FOUND="T"
        fi
        kold=$k
      done
    done


   echo i=$i
   exit 1

  fi
done
rm -f .dir* .tmp .file1 .file2
