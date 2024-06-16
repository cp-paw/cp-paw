#!/bin/bash
#--  help message --------------------------------------------------------------
export USAGE=""
USAGE="${USAGE}\n usage: paw_checkpoint [create|recover] name"
USAGE="${USAGE}\n"
USAGE="${USAGE}\n create:\t copy files from current directory to checkpoint"
USAGE="${USAGE}\n recover:\t delete all files in the current directory and"
USAGE="${USAGE}\n \t\t copy content from checkpoint into current directory"
USAGE="${USAGE}\n name\t\t name of the checkpoint directory"

#-------------------------------------------------------------------------------
#--  check input syntax 
#-------------------------------------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "error in $0: Illegal number of parameters"
  echo -e ${USAGE}
  exit 1
fi

if [[ $1 != create && $1 != recover ]] ; then
  echo "error in $0: first argument must be \"create\" or \"recover\""
  echo "first argument is $1"
  echo -e ${USAGE}
  exit 1
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
OP="$1"
DIRNAME=$2

case "$OP" in
  create)
    #check if directory $DIRNAME already exists
    if [[ -d $DIRNAME ]] ; then
      echo "WARNING: the directory $DIRNAME already exists, doing NOTHING!" 
      exit
    else
      mkdir "$DIRNAME" # create directory
      /bin/cp * "$DIRNAME" #copy all files
    fi
  ;;
  recover)
    #check if directory $DIRNAME exists
    if [ -d "$2" ] ;  then
      #ask user if he really wants to
      echo "Do you really want to delete all files (not directories) in the current directory and copy the files from $DIRNAME into it? (y/n)"
      read answer
      if [ "$answer" = "y" ];  then
        #remove all files in the current directory
        /bin/rm *
        #copy all files from $DIRNAME
        /bin/cp $DIRNAME/* .
        echo "You now have the all the files from the directory $DIRNAME in the current directory."
      else
        echo "DOING NOTHING."
        exit
      fi
    else
      echo "WARNING: directory $DIRNAME does NOT exist, doing NOTHING!" 
      exit
    fi
  ;;
esac


