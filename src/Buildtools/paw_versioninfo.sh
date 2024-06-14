#!/bin/bash
export THISDIR=$(pwd)
###########################################################################
#
#  Purpose:
#     collect information about the last commit from git repository
#     and write information to standard out
#
#     if git is not installed or command is executed outside a working tree
#     the error code $? is set to 1
#     and /dev/stout remains empty  (except for verbose mode)
#
###########################################################################
export OUT=/dev/stdout   #redirect to standard out
export VERBOSE=F

#----------------------------------------------------
# default values
#----------------------------------------------------
export REMOTE='unknown'
export BRANCH='unknown'
export HASH='unknown'
export SHORTREVISIONNUMBER='unknown'
export AUTHOR='unknown'
export DATE='unknown'
#
#----------------------------------------------------
# check if git is installed. string is empty if not
#----------------------------------------------------
GITVERSION=$(git version 2>/dev/null)
#
if [[ $VERBOSE = T ]] ; then
  if [[ -z $GITVERSION ]] ; then
    echo 'git is not installed'                                              >&2
  else
    echo 'git is installed in version ' $GITVERSION
  fi
fi
if [[ -z ${GITVERSION} ]] ; then exit 1; fi
#
#----------------------------------------------------
# check if inside a git working tree 
#----------------------------------------------------
if [[ -n $(git rev-parse --git-dir 2>/dev/null)  ]] ; then 
  #----------------------------------------------------
  # location of the remote repository
  #----------------------------------------------------
  # REMOTE=$(git remote -v | awk 'BEGIN { FS = " " } ; { print $2 }' | uniq | tr '\n' ';')
  # The command above lists all remotes, while the following lists only the first
  #
  REMOTE=$(git remote -v | awk 'BEGIN { FS = " " } ; { print $2 }' | head -n 1)
  if [[ $VERBOSE = T ]] ; then
    echo REMOTE=$REMOTE
  fi
  #----------------------------------------------------
  # identify active branch
  #----------------------------------------------------
  BRANCH=$(git branch | grep "*" | sed "s/* //")
  if [[ $VERBOSE = T ]] ; then
    echo BRANCH=$BRANCH
  fi
  #----------------------------------------------------
  # short revisionnumber
  #----------------------------------------------------
  export SHORTREVISIONNUMBER=$(git rev-list $BRANCH | wc -l)
  if [[ $VERBOSE = T ]] ; then
    echo SHORTREVISIONNUMBER=$SHORTREVISIONNUMBER
  fi
  #----------------------------------------------------
  # hash of the last commit
  #----------------------------------------------------
  export HASH=$(git log -1 | head -n 1| cut -d " " -f2)
  if [[ $VERBOSE = T ]] ; then
    echo HASH=$HASH
  fi
  #----------------------------------------------------
  # number of changes since last commit
  #----------------------------------------------------
  export NUMCHANGES=$(git diff |grep @@ |wc -l)
  if [[ $VERBOSE = T ]] ; then
    echo NUMCHANGES=$NUMCHANGES
  fi
  #----------------------------------------------------
  # author
  #----------------------------------------------------
  AUTHOR=$(git log -1 | grep "Author: " | head -n 1 | sed "s/Author: //g")
  if [[ $VERBOSE = T ]] ; then
    echo AUTHOR=$AUTHOR
  fi
  #----------------------------------------------------
  # date
  #----------------------------------------------------
  DATE=$(git log -1 | grep "Date: " | head -n 1 | sed "s/Date:   //g") 
  if [[ $VERBOSE = T ]] ; then
    echo DATE=$DATE
  fi
else
  echo "error in $0: not under git control"                                  >&2
  exit 2   # not under git control
fi
#---------------------------------------------------------------------------
# write to version info to $OUT
#---------------------------------------------------------------------------
echo "REMOTE=       '$REMOTE'"                    >  $OUT
echo "BRANCH=       '$BRANCH'"                    >> $OUT
echo "HASH=         '$HASH'"                      >> $OUT
echo "SHORTREVISIONNUMBER='$SHORTREVISIONNUMBER'" >> $OUT
echo "AUTHOR=       '$AUTHOR'"                    >> $OUT
echo "COMMITDATE=   '$DATE'"                      >> $OUT 
echo "NUMCHANGES=    $NUMCHANGES"                 >> $OUT 
echo "COMPILEDATE=  '$(date)'"                    >> $OUT 
echo "COMPILEPERSON='$(whoami)'"                  >> $OUT 
