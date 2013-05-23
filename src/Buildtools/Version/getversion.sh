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
    ##all remotes--> too long
    VER1=`git remote -v | awk 'BEGIN { FS = " " } ; { print $2 }' | uniq | tr '\n' ';'`
    ##only remote that hast tracking branch
    #branch=$(git symbolic-ref HEAD)
    #branch=${branch##refs/heads/}
    #remote=$(git config "branch.${branch}.remote")
    #remoteBranch=$(git config "branch.${branch}.merge")
    #remoteBranch=${remoteBranch##refs/heads/}
    #VER1=`git remote show -n ${remote:?} | grep "Fetch URL:" | sed "s/[ ]*Fetch URL://g" | tr -d '\n'`
    #banchname
    VER2=`git branch | grep "*" | sed "s/* //"`
    #short revision number
    VER3=`git rev-list $VER2 | wc -l`
    #hash of last commit
    #VER4=`git log -1 --format="%H"` # --> not working on git 1.5.5.1
    VER4=`git log -1 | head -n 1| cut -d " " -f2`
    #number of files changed since last commit
    #VER5=`git status --porcelain ../../../../src | wc -l` # --> not working on git 1.5.5.1
    VER5=`git diff ../../../../src | wc -l`
    #AUTHOR=`git log -1 --format="%an (%ae)"` # --> not working on git 1.5.5.1
    AUTHOR=`git log -1 | grep "Author: " | head -n 1 | sed "s/Author: //g"`
    #DATE=`git log -1 --format="%aD"` #--> not working on git 1.5.5.1
    DATE=`git log -1 | grep "Date: " | head -n 1 | sed "s/Date:   //g"`
    VERINF="branch: $VER2; remote: $VER1"
    VERREV="revision $VER3; $VER4; $VER5 changes since last commit"
    VERAUT="last commit by $AUTHOR; compiled by `whoami` on `hostname`"
    VERDAT="last commit at $DATE; compiled at `date`"    
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
VERINF2=`echo "$VERINF" | sed 's/\//\\\\\//g'| cut  -c -250`
VERREV2=`echo "$VERREV" | sed 's/\//\\\\\//g'| cut  -c -250`
VERAUT2=`echo "$VERAUT" | sed 's/\//\\\\\//g'| cut  -c -250`
VERDAT2=`echo "$VERDAT" | sed 's/\//\\\\\//g'| cut  -c -250`

#modify in file
sed -i "s/_VERINF/${VERINF2}/g" $PAWVERSIONFILE
sed -i "s/_VERREV/${VERREV2}/g" $PAWVERSIONFILE
sed -i "s/_VERAUT/${VERAUT2}/g" $PAWVERSIONFILE
sed -i "s/_VERDAT/${VERDAT2}/g" $PAWVERSIONFILE


