#!/bin/bash
export THISDIR=$(pwd)
#############################################################################
#
#   make_release_package options
#
#   purpose:
#
#      create a distribution package from the current working copy
#
#   options:
#      -i version identifier (e.g. v2024.1)
#      -b name of the resulting tarball, optional (e.g. cppaw_v2024.1.tgz) 
#      -p parameter file (e.g. parms.blub.bla)
#      -h help information
#
#    this script has to be called in the root of the working copy, 
#    i.e. in the directory containing src, Docs,...
#
#   example: 
#     sh src/Buildtools/Make_distribution_package/make_release_package.sh \
#        -i v2024.1 -p parms.linux
#
# all the cool stuff written by Robert Schade
# modified by Peter Bloechl
#############################################################################
# NO MORE: multiple parms-files are possible, the fist one is used to create docs, so it should be working on the current machine

#as temporary storage $WD is used 

#-------------------------------------------------------------------------------
# help message
#-------------------------------------------------------------------------------
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \t src/Buildtools/Make_distribution_package/make_release_package.sh tarball versionid parmfile\n"
USAGE="$USAGE \n"
USAGE="$USAGE Purpose:\n"
USAGE="$USAGE \t produces a tarball of the cp-paw distribution \
                 with builtin version information\n"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
# USAGE="$USAGE \t tarball \t name of the tar-ball with the distribution \
#                             e.g. CP-PAW_V2023.1.tar.gz\n"
# USAGE="$USAGE \t versionid \t version identifier, e.g. v2023.1\n"
# USAGE="$USAGE \t parmfile \t name of the parmfile for the\
#                  configuration of the CP-PAW\n"
USAGE="$USAGE \t -b \t name of the tarball of the distribution with full path\n"
USAGE="$USAGE \t    \t without path it is relative to /tmp\n"
USAGE="$USAGE \t -i \t optional version identifier, e.g. v2023.1\n"
USAGE="$USAGE \t -p \t parmfile\n"
USAGE="$USAGE \t -h \t print this help message \n"
# USAGE="$USAGE \t -0: dry-run (creates files but does not run jobs)\n"
# USAGE="$USAGE \t -v: verbose\n"
USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
# resolve command line arguments
#-------------------------------------------------------------------------------
export TYPE='RELEASE'
#export TYPE='DEVELOPMENT'
export TARBALL=""
export VERSIONID=""
export PARMFILE=""

OPTSTRING=":hv0b:i:p:"
OPTIND=0
while getopts "${OPTSTRING}" OPT  ; do
  case $OPT in
    b) TARBALL="$OPTARG" #name of the tarball to be created
      ;;
    i) VERSIONID="$OPTARG" # version id
      ;;
    p) PARMFILE="$OPTARG" # name of the parmfile
      ;;
    0)   #nothing:
      DRYRUN=yes
#      set -n
      ;;
    v)   #verbose
      VERBOSE=yes
#      set -v
#      set -x
      ;;
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
# report command line arguments and catch missing input
#-------------------------------------------------------------------------------
#  default values for tarball
if [[ -z $TARBALL && -n $VERSIONID ]] ; then
  TARBALL="${THISDIR}/CP-PAW_"${VERSIONID}".tar.gz"
fi

# default value for parmfile 
if [[ -z $PARMFILE ]] ; then
  if [[ -f parms.in_use ]] ; then
     PARMFILE="parms.in_use"
  fi
fi

# report current seting
echo -e "name of the tarball:  $TARBALL"
echo -e "version id.........:  $VERSIONID"
echo -e "parmfile to be used:  $PARMFILE"

# error message for missing input
if [[ -z $TARBALL || -z $VERSIONID || -z $PARMFILE ]] ; then 
  if [[ -z $TARBALL ]] ; then 
    echo 'error: name of tarball (-b) not specified'
  fi
  if [[ -z $VERSIONID ]] ; then 
    echo 'error: version id (-i) not specified'
  fi
  if [[ -z $PARMFILE ]] ; then 
    echo 'error: parmfile (-p) not specified'
  fi
  echo 'error'
  exit 1
fi

# error exit for missing parmfile
if [[ ! -f  $PARMFILE ]] ; then
  echo 'error: parmfile $PARMFILE does not exist'
  exit 1
fi

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
if [[ -f ${TARBALL} ]] ; then 
  rm ${TARBALL}
fi

export WORKDIR="/tmp/cp-paw_$VERSIONID"
if [[ -d $WORKDIR ]] ; then rm -rf $WORKDIR; fi
mkdir -p $WORKDIR
cp -r * $WORKDIR
cp -r .git $WORKDIR
cp -v $PARMFILE $WORKDIR
cd $WORKDIR
#
#-------------------------------------------------------------------------------
#  construct documentation and clean $WORKDIR
#-------------------------------------------------------------------------------
echo "configuring cppaw distribution....."
./configure --with-parmfile=$(basename $PARMFILE)
# construct cppaw_version.info in the top directory 
echo "making cppaw_version.info....."
make version 1>/dev/null 2>&1
if [[ $? -ne 0 ]] ; then echo "error: no version information" ; exit 1 ; fi
if [[ -s cppaw_version.info ]] ; then
  echo "RELEASE= '$VERSIONID'" >> cppaw_version.info
fi
echo "making documentation....."
# construct documentation in doc directy
make docs 1>/dev/null 2>&1
if [[ $? -ne 0 ]] ; then echo "latex compilation error" ; exit 1 ; fi
make clean
rm -rf .git

#-------------------------------------------------------------------------------
# make a tarball of $WORKDIR and remove $WORKDIR
#-------------------------------------------------------------------------------
cd ..
tar -cvzf $TARBALL $(basename $WORKDIR)
rm -rf $WORKDIR
