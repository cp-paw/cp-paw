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
#      -h help information
#
#    this script has to be called in the root of the working copy, 
#    i.e. in the directory containing src, Docs,...
#
#   example: 
#     sh src/Buildtools/Make_distribution_package/make_release_package.sh \
#        -i v2024.1 
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
USAGE="$USAGE \t -b \t name of the tarball of the distribution with full path\n"
USAGE="$USAGE \t    \t without path it is relative to /tmp\n"
USAGE="$USAGE \t -i \t optional version identifier, e.g. v2023.1\n"
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

OPTSTRING=":hv0b:i:p:"
OPTIND=0
while getopts "${OPTSTRING}" OPT  ; do
  case $OPT in
    b) TARBALL="$OPTARG" #name of the tarball to be created
      ;;
    i) VERSIONID="$OPTARG" # version id
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

# report current seting
echo -e "name of the tarball:  $TARBALL"
echo -e "version id.........:  $VERSIONID"

# error message for missing input
if [[ -z $TARBALL || -z $VERSIONID ]] ; then 
  if [[ -z $TARBALL ]] ; then 
    echo 'error: name of tarball (-b) not specified'
  fi
  if [[ -z $VERSIONID ]] ; then 
    echo 'error: version id (-i) not specified'
  fi
  echo 'error'
  exit 1
fi

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
if [[ -f ${TARBALL} ]] ; then 
  rm ${TARBALL}
fi
#
#-------------------------------------------------------------------------------
#  construct version information. will be used during installation
#-------------------------------------------------------------------------------
# construct cppaw_version.info in the top directory 
echo "collect version information....."
echo "RELEASE= '$VERSIONID'" > cppaw_version.info
./src/Buildtools/paw_versioninfo.sh >> cppaw_version.info
#
#-------------------------------------------------------------------------------
#  create temporary working directory for the construction of the tar ball
#-------------------------------------------------------------------------------
export TMPDIR=$(mktemp -d)
export WORKDIR=${TMPDIR}/cppaw
if [[ ! -d ${WORKDIR} ]] ; then mkdir ${WORKDIR}; fi
cp -r * $WORKDIR
#cp -r .git $WORKDIR
cd $WORKDIR
#
#-------------------------------------------------------------------------------
#  construct documentation and clean $WORKDIR
#-------------------------------------------------------------------------------

# echo "ensure to use the local directory"
# export PATH=${THISDIR}/src/Buildtools/:${THISDIR}/src/Tools/Scripts/:$PATH

echo "install cppaw distribution....."
./paw_install
if [[ $? -ne 0 ]] ; then 
  echo "paw_installation failed" 
  echo "installation in workdirectory ${WORKDIR}" 
  ls -al ./
  exit 1 
fi

echo "collect parms_in_use....."
./bin/fast/paw_fast.x -p 1>/dev/null

#-------------------------------------------------------------------------------
# cleanup
#-------------------------------------------------------------------------------
rm -rf ${WORKDIR}/bin
rm -rf ${WORKDIR}/.git
rm -rf ${WORKDIR}/err_*
rm -rf ${WORKDIR}/out_*

#-------------------------------------------------------------------------------
# make a tarball of $WORKDIR and remove $WORKDIR
#-------------------------------------------------------------------------------
cd ..
tar -cvzf $TARBALL $(basename $WORKDIR)
rm -rf $WORKDIR


