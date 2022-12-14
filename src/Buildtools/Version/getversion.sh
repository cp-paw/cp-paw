#!/bin/bash
PAWVERSIONFILE="$1"
echo "==================getversion.sh start==================================="
echo "File to modify: $PAWVERSIONFILE"
echo "Current directory=$(pwd)"

srcdir="../../../../src"
if [ "$PAWVERSIONFILE" = "" ] ;
then
  srcdir="src"
fi

VERTYP="'UNKNOWN VERSION'"
VERINF="'unknown'"
VERREV="'unknown'"
VERAUT="'unknown'"
VERDAT="'unknown'"


#check if version info file exists
verinfo="$srcdir/version.info"
if [ -e "${verinfo}" ]
then
  VERTYP=`sed -n '1p' ${verinfo}`
  if [ "$VERTYP" = "'DEVELOPMENT VERSION'" ];
  then
    VER1=`sed -n '2p' ${verinfo}`
    VER2=`sed -n '3p' ${verinfo}`
    VER3=`sed -n '4p' ${verinfo}`
    VER4=`sed -n '5p' ${verinfo}`
    VER5=`sed -n '6p' ${verinfo}`
    AUTHOR=`sed -n '7p' ${verinfo}`
    DATE=`sed -n '8p' ${verinfo}`

    VERINF="'branch: $VER2; remotes: $VER1'"
    VERREV="'revision $VER3; $VER4; $VER5 changes since last commit'"
    VERAUT="'last commit by $AUTHOR at $DATE'"
    VERDAT="'compiled by `whoami` on `hostname` at `date`'"    
  elif [ "$VERTYP" = "'RELEASE VERSION'" ];
  then
    VERSION=`sed -n '2p' ${verinfo}`
    DATE=`sed -n '3p' ${verinfo}`
    AUTHOR=`sed -n '4p' ${verinfo}`

    VERINF="'$VERSION'"
    VERREV="'$DATE'"
    VERAUT="'$AUTHOR'"
    VERDAT="'compiled by `whoami` on `hostname` at `date`'"    
  fi
else
  #check if in git installed
  git --version 2>&1 >/dev/null # improvement by tripleee
  GIT_IS_AVAILABLE=$?
  if [ "$GIT_IS_AVAILABLE" = 0 ]; then
    echo "git installed: YES"    
    #git ist installed
    if git rev-parse --git-dir > /dev/null 2>&1; then
      #git is installed and we are in a git repo
      echo "inside git repo: YES"    
      ##all remotes--> too long
      # VER1 is a ;-separated list of the remote repositories
      VER1=`git remote -v | awk 'BEGIN { FS = " " } ; { print $2 }' | uniq | tr '\n' ';'`
      ##only remote that hast tracking branch
      #branch=$(git symbolic-ref HEAD)
      #branch=${branch##refs/heads/}
      #remote=$(git config "branch.${branch}.remote")
      #remoteBranch=$(git config "branch.${branch}.merge")
      #remoteBranch=${remoteBranch##refs/heads/}
      #VER1=`git remote show -n ${remote:?} | grep "Fetch URL:" | sed "s/[ ]*Fetch URL://g" | tr -d '\n'`
      #banchname
      # VER2 is the name of the actual branch (with an unintended line number in front)
      VER2=`git branch | grep "*" | sed "s/* //"`
      #short revision number
      VER3=`git rev-list $VER2 | wc -l`
      #hash of last commit
      #VER4=`git log -1 --format="%H"` # --> not working on git 1.5.5.1
      VER4=`git log -1 | head -n 1| cut -d " " -f2`
      #number of files changed since last commit
      #VER5=`git status --porcelain ../../../../src | wc -l` # --> not working on git 1.5.5.1
      VER5=`git diff $srcdir | wc -l`
      #AUTHOR=`git log -1 --format="%an (%ae)"` # --> not working on git 1.5.5.1
      AUTHOR=`git log -1 | grep "Author: " | head -n 1 | sed "s/Author: //g"`
      #DATE=`git log -1 --format="%aD"` #--> not working on git 1.5.5.1
      DATE=`git log -1 | grep "Date: " | head -n 1 | sed "s/Date:   //g"`
      VERTYP="'DEVELOPMENT VERSION'"
      VERINF="'branch: $VER2; remotes: $VER1'"
      VERREV="'revision $VER3; $VER4; $VER5 changes since last commit'"
      VERAUT="'last commit by $AUTHOR at $DATE'"
      VERDAT="'compiled by `whoami` on `hostname` at `date`'"    
    else
      #git is installed, but we are not in a git repo
      echo "inside git repo: NO"    
      VERTYP="'UNKNOWN VERSION'"
      VERINF="'unknown'"
      VERREV="'unknown'"
      VERAUT="'unknown'"
      VERDAT="'compiled by `whoami` on `hostname` at `date`'"    
    fi
  else
    #git is not installed
    echo "git installed: NO"    
    VERTYP="'UNKNOWN VERSION'"
    VERINF="'unknown'"
    VERREV="'unknown'"
    VERAUT="'unknown'"
    VERDAT="'compiled by `whoami` on `hostname` at `date`'"    
  fi
fi

#echo "VERINF=$VERINF2"
#echo "VERREV=$VERREV2"
#echo "VERAUT=$VERAUT2"
#echo "VERDAT=$VERDAT2"

#only change module file if $PAWVERSIONFILE is not empty
if [ "$PAWVERSIONFILE" = "" ]
then
  echo "not changing version information in paw-code, just writing it to $verinfo"
  echo "pwd=$(pwd)"
  echo "verinfo=${verinfo}"
  echo "'DEVELOPMENT VERSION'" > ${verinfo}
  echo "$VER1" >> ${verinfo}
  echo "$VER2" >> ${verinfo}
  echo "$VER3" >> ${verinfo}
  echo "$VER4" >> ${verinfo}
  echo "$VER5" >> ${verinfo}
  echo "$AUTHOR" >> ${verinfo}
  echo "$DATE" >> ${verinfo}
else

  #fixing line end if length is multiple of separation length
  l=`echo "$VERTYP" | wc -c`
  m=`echo "($l-1)%50==0" | bc`
  if [ "$m" = "1" ];
  then
    VERTYP=`echo "$VERTYP" | sed "s/'$/ '/g"`
  fi
  l=`echo "$VERINF" | wc -c`
  m=`echo "($l-1)%50==0" | bc`
  if [ "$m" = "1" ];
  then
    VERINF=`echo "$VERINF" | sed "s/'$/ '/g"`
  fi
  l=`echo "$VERREV" | wc -c`
  m=`echo "($l-1)%50==0" | bc`
  if [ "$m" = "1" ];
  then
    VERREV=`echo "$VERREV" | sed "s/'$/ '/g"`
  fi
  l=`echo "$VERAUT" | wc -c`
  m=`echo "($l-1)%50==0" | bc`
  if [ "$m" = "1" ];
  then
    VERAUT=`echo "$VERAUT" | sed "s/'$/ '/g"`
  fi
  l=`echo "$VERDAT" | wc -c`
  m=`echo "($l-1)%50==0" | bc`
  if [ "$m" = "1" ];
  then
    VERDAT=`echo "$VERDAT" | sed "s/'$/ '/g"`
  fi

  # $(command) is the output of teh command as text string
  # \n is a newline character to brak over-long lines in the fortran code
  # // concatenates the strings on the two separate lines
  VERTYP2=$(echo "$VERTYP" |  sed "s|.\{50\}|&\' \/\/ \& \\\n  \'|g")
  VERINF2=$(echo "$VERINF" |  sed "s|.\{50\}|&\' \/\/ \& \\\n  \'|g")
  VERREV2=$(echo "$VERREV" |  sed "s|.\{50\}|&\' \/\/ \& \\\n  \'|g")
  VERAUT2=$(echo "$VERAUT" |  sed "s|.\{50\}|&\' \/\/ \& \\\n  \'|g")
  VERDAT2=$(echo "$VERDAT" |  sed "s|.\{50\}|&\' \/\/ \& \\\n  \'|g")
  echo "MODULE VERSION_MODULE" > $PAWVERSIONFILE.tmp
  echo -e "CHARACTER(256):: VERTYP=$VERTYP2" >> $PAWVERSIONFILE.tmp
  echo -e "CHARACTER(256):: VERINF=$VERINF2" >> $PAWVERSIONFILE.tmp
  echo -e "CHARACTER(256):: VERREV=$VERREV2" >> $PAWVERSIONFILE.tmp
  echo -e "CHARACTER(256):: VERAUT=$VERAUT2" >> $PAWVERSIONFILE.tmp
  echo -e "CHARACTER(256):: VERDAT=$VERDAT2" >> $PAWVERSIONFILE.tmp
  echo "END MODULE VERSION_MODULE" >> $PAWVERSIONFILE.tmp
  cat $PAWVERSIONFILE.tmp $PAWVERSIONFILE > $PAWVERSIONFILE.tmp2
  rm $PAWVERSIONFILE.tmp
  mv $PAWVERSIONFILE.tmp2 $PAWVERSIONFILE
fi
echo "==================getversion.sh end==================================="


