#!/bin/bash
#
#  paw_getsrc.sh pawexecutable pawsrc.tgz
#
#  pawexecutable is for example paw_fast.x
#  pawsrc.tgz is the name of the zipped tar file containing the source files
#  
#  inspect with: tar tvzf pawsrc.tgz
#  unpack with:  tar xvzf pawsrc.tgz
#  
# this script extracts the embedded source code from a cppaw binary
# which is in macho-format if you have a binary that has been created
# on linux/BSD please use getsrc_elf.sh call with: bash getsrc_elf.sh
# $PAWBINARY $DESTINATION PAWBINARY is for example paw_fast.x
# DESTINATION is the location for the output of the source-tar-archive
#
export THIS=$(basename "$0")
#-------------------------------------------------------------------------------
# help message
#-------------------------------------------------------------------------------
export USAGE="Usage of $0 \n"
USAGE="$USAGE \n"
USAGE="$USAGE \t${THIS}\n"
USAGE="$USAGE \n"
USAGE="$USAGE Purpose:\n"
USAGE="$USAGE \t extract embedded sourceblob from a paw executable"
USAGE="$USAGE \n"
USAGE="$USAGE Options:\n"
USAGE="$USAGE \t -h: print this help message \n"
USAGE="$USAGE \t -x: name of the paw_executable \n"
USAGE="$USAGE \t -o: name of the sourceblob \n"
USAGE="$USAGE \t -0: dry-run (creates files but does not run jobs)\n"
USAGE="$USAGE \t -v: verbose\n"
USAGE="$USAGE \n"
USAGE="$USAGE Example:\n"
USAGE="$USAGE \t ${THIS} -x paw_fast.x -o paw_src.tgz \n"
USAGE="$USAGE \n"
#-------------------------------------------------------------------------------
#  resolve argument list
#-------------------------------------------------------------------------------
export PAWX=$(which paw_fast.x)
export SRCBLOB=pawsrc.tgz
#
OPTSTRING=":hv0x:o:"
OPTIND=0
while getopts "${OPTSTRING}" OPT  ; do
  case $OPT in
    x) PAWX=$OPTARG #set executable. Default ist $(which paw_fast.x)
      ;;
    o) SRCBLOB="$OPTARG" #set name of the source blob. (default is pawsrc.tgz)
      ;;
    0)   #nothing:
      set -n
      ;;
    v)   #verbose
      set -v
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
#
# avoid overwriting an existing file
if [[ ! -f $PAWX ]] ; then 
  echo "error: executable $PAWX does not exist" >&2
  echo "error in $0" >&2
  exit 1
fi
if [[ -e $SRCBLOB ]] ; then 
  echo "error: target file $SRCBLOB for the source blob already exists" >&2
  echo "error in $0" >&2
  exit 1
fi
#
#===================================================================
#  make two temporary files
#===================================================================
# check if environment variable TMPDIR is set
if [[ -z ${TMPDIR} ]] ; then TMPDIR=/tmp ; fi
if [[ ! -d ${TMPDIR} ]] ; then
  echo "error in $0: temp directory does not exist" >&2
  echo "TMPDIR=$TMPDIR" >&2
  exit 1
fi
TMPFILE1=$(mktemp ${TMPDIR}/pawsrcblobtmp.XXXXXX)
RC=$?
if [[ $RC -ne 0 ]]; then
  echo "error in $0: Can't create temp file, exiting..."
  exit 1
fi
TMPFILE2=$(mktemp ${TMPDIR}/pawsrcblobtmp.XXXXXX)
RC=$?
if [[ $RC -ne 0 ]]; then
  echo "error in $0: Can't create temp file, exiting..."
  exit 1
fi
#
#===================================================================
#  make two temporary files
#===================================================================
case $(uname) in
  Darwin) # Apple environment
    # otool extracts a section named paw_srcblob_bin from the paw executable.
    # this section contains the tar file with the source code.
    # the tar command creating the source code section is
    #       tar -czf paw_srcblob ${SRCDIR} ${THISDIR}/parms.in_use
    # the otool produces a hexdump of this source code section
    # 
    otool $PAWX -s binary pawsrcblob_bin > $TMPFILE1
    #
    # removes irrelevant information (first two columns) and combines the
    # rest into a single hexdump file
    #
    tail +3 $TMPFILE1 | awk 'BEGIN { FS = " " } ; { print $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 }' > $TMPFILE2
    #
    # xxd converts a hexdump  file into a binary
    # 
    if [[ $(uname -p) -eq 'arm' ]] ; then # for processor type (ARM)
      #
      # Converts little endian in big endian or vice versa
      #
      xxd -r -p $TMPFILE2 $TMPFILE1
      hexdump -v -e '1/4 "%08x"'  $TMPFILE1 >$TMPFILE2
    fi 
    xxd -r -p ${TMPFILE2} ${SRCBLOB}
    ;;
  Linux) 
    start=$(objdump -t $PAWX \
             | grep _binary_paw_srcblob_start \
             | awk 'BEGIN { FS = " " } ; { print $1 }' \
             | sed "s/^0*//g" \
           )
    end=$(objdump -t $PAWX \
             | grep _binary_paw_srcblob_end \
             | awk 'BEGIN { FS = " " } ; { print $1 }' | sed "s/^0*//g" \
           )
    size=$(objdump -t $PAWX \
             | grep _binary_paw_srcblob_size \
             | awk 'BEGIN { FS = " " } ; { print $1 }' | sed "s/^0*//g" \
           )
    startdata=$(objdump -h $PAWX \
                 | grep " .data " \
                 | awk 'BEGIN { FS = " " } ; { print $4 }' | sed "s/^0*//g"
               )
    startdatareal=$(objdump -h $PAWX \
                    | grep " .data " \
                    | awk 'BEGIN { FS = " " } ; { print $6 }' | sed "s/^0*//g" \
                   )
    startblobreal=$(echo "$[0x$start]-$[0x$startdata]+$[0x$startdatareal]" | bc)
    dd if=$PAWX of=$SRCBLOB bs=1 count=$[0x$size] skip=$startblobreal
    ;;
  *)
    echo "error in $0: unknown operating system" 
    echo "Known are 'Darwin' and 'Linux' " 
    exit 1
    ;;
esac
#
# clean up
#
rm $TMPFILE1
rm $TMPFILE2
