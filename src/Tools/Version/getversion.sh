#!/bin/bash
PAWVERSIONFILE="$1"
echo "File to modify: $PAWVERSIONFILE"

VERINF=""
VERREV=""
VERAUT=""
VERDAT=""

#check if in git installed
git --version 2>&1 >/dev/null # improvement by tripleee
GIT_IS_AVAILABLE=$?
if [ "$GIT_IS_AVAILABLE" = 0 ]; then
  echo "git installed: YES"    
  #git ist installed
  if git rev-parse --git-dir > /dev/null 2>&1; then
    #git is installed and we are in a git repo
    echo "inside git repo: YES"    
    #remotes
    VER1=`git remote -v | awk 'BEGIN { FS = " " } ; { print $2 }' | uniq | tr '\n' ';'`
    #banchname
    VER2=`git branch | grep "*" | sed "s/* //"`
    #short revision number
    VER3=`git rev-list $VER2 | wc -l`
    #hash of last commit
    VER4=`git log -1 --format="%H"` 
    #number of files changed since last commit
    VER5=`git status --porcelain ../../../../src | wc -l` 
    VERINF="remotes: $VER1 branch: $VER2"
    VERREV="r$VER3; $VER4; $VER5 changes since last commit"
    VERAUT="last commit by `git log -1 --format="%an (%ae)"`; compiled by `whoami` on `hostname`"
    VERDAT="last commit at `git log -1 --format="%aD"`; compiled at `date`"    
  else
    #git is installed, but we are not in a git repo
    echo "inside git repo: NO"    
    VERINF="unknown"
    VERREV="unknown"
    VERAUT="compiled by `whoami` on `hostname`"
    VERDAT="compiled at `date`"
  fi
else
  #git is not installed
  echo "git installed: NO"    
  VERINF="unknown"
  VERREV="unknown"
  VERAUT="compiled by `whoami` on `hostname`"
  VERDAT=`date`
fi

#echo "VERINF=$VERINF"
#echo "VERREV=$VERREV"
#echo "VERAUT=$VERAUT"
#echo "VERDAT=$VERDAT"
VERINF2=`echo "$VERINF" | sed 's/\//\\\\\//g'`
VERREV2=`echo "$VERREV" | sed 's/\//\\\\\//g'`
VERAUT2=`echo "$VERAUT" | sed 's/\//\\\\\//g'`
VERDAT2=`echo "$VERDAT" | sed 's/\//\\\\\//g'`

#modify in file
sed -i "s/_VERINF/${VERINF2}/g" $PAWVERSIONFILE
sed -i "s/_VERREV/${VERREV2}/g" $PAWVERSIONFILE
sed -i "s/_VERAUT/${VERAUT2}/g" $PAWVERSIONFILE
sed -i "s/_VERDAT/${VERDAT2}/g" $PAWVERSIONFILE



