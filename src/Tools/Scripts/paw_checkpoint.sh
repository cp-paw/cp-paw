#!/bin/bash
#--  help message --------------------------------------------------------------
export USAGE=""
USAGE="${USAGE}\n usage: paw_checkpoint option"
USAGE="${USAGE}\n"
USAGE="${USAGE}\n\t -c name:\t create snapshop directory name"
USAGE="${USAGE}\n\t\t\t 1) create a directory name"
USAGE="${USAGE}\n\t\t\t 2) copy files (no directories) from current "
USAGE="${USAGE} to snapshot directory"
USAGE="${USAGE}\n\t -r name:\t recover snapshot name"
USAGE="${USAGE}\n\t\t\t 1) erase all files (no directories) in current directory"
USAGE="${USAGE}\n\t\t\t 2) copy files from snapshot into current directory"
USAGE="${USAGE}\n"

#-------------------------------------------------------------------------------
#--  resolve options
#-------------------------------------------------------------------------------
export DRYRUN="false"
export VERBOSE="false"
export OP=""
export DIRNAME=""
OPTSTRING=":hc:r:"
while getopts "${OPTSTRING}" OPT  ; do
  case $OPT in
    c) DIRNAME="$OPTARG" 
       if [[ -n $OP ]] ; then 
         echo "error in $0: operation can only be defined once"
         echo "OPT=$OPT and operation is $OP"
         exit 1
       fi
       OP="create"
       ;;
    r) DIRNAME="$OPTARG" 
       if [[ -n $OP ]] ; then 
         echo "error in $0: operation can only be defined once"
         echo "OPT=$OPT and operation is $OP"
         exit 1
       fi
       OP="recover"
       ;;
#     0)   #nothing:
#        DRYRUN="true"
# #      set -n
#        ;;
    # v)   #verbose
    #    VERBOSE="true"
#       set -v
#      set -x
#      ;;
    h)   # help
      echo -e $USAGE
      exit 0
      ;;
    \?)   # unknown option (placed into OPTARG, if OPTSTRING starts with :)
      echo "error in $0" >&2
      echo "invalid option -$OPTARG" >&2
      echo "retrieve argument list with:" >&2
      echo "$0 -h" >&2
      exit 1
      ;;
    :)    # no argument passed to option requiring one
      echo "error in $0" >&2
      echo "option -$OPTARG requires an additional argument" >&2
      exit 1
      ;;  
  esac
done

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
case "$OP" in
  create)
    #check if directory $DIRNAME already exists
    if [[ -d ${DIRNAME} ]] ; then
      echo "WARNING: the directory $DIRNAME already exists, doing NOTHING!" 
      exit 1
    else
      if [[ $DRYRUN != true ]] ; then
        mkdir "$DIRNAME" # create directory
        /bin/cp * "$DIRNAME" #copy all files
      fi
    fi
    ;;
#
  recover)
    #check if directory $DIRNAME exists
    if [ -d "${DIRNAME}" ] ;  then
      #ask user if he really wants to
      TEXT=""
      TEXT="${TEXT} Do you really want to delete all files (not directories)"
      TEXT="${TEXT} in the current directory\n"
      TEXT="${TEXT} and copy the files from $DIRNAME into it? (y/n)"
      echo -e "$TEXT"
      read answer
      if [ "$answer" = "y" ];  then
        #remove all files in the current directory
        if [[ $DRYRUN != true ]]; then
          /bin/rm *
          #copy all files from $DIRNAME
          /bin/cp $DIRNAME/* .
          echo "Now you have the all files from the directory $DIRNAME"\
             " in the current directory."
        fi
      else
        echo "DOING NOTHING."
        exit 1
      fi
    else
      echo "WARNING: directory $DIRNAME does NOT exist, doing NOTHING!" 
      exit
    fi
  ;;
#
  *) #do nothing
    ;; 
esac


