#!/bin/sh
#*******************************************************************
#*******************************************************************
#****                                                           ****
#****  creates a mpeg movie file from a set of gif files        ****
#****                                                           ****
#****  $1 is the root of the gif files which have names         ****
#****      $1.num.gif where num is a number starting at         ****
#****      include path in $1 such as ./movie                   ****
#****                                                           ****
#****  it uses the mpeg_encode executable defined as            ****
MPEG_ENCODE=/u/blo/bin/mpeg_encode 
GIFTOPPM=/u/blo.2/Resources/giftoppm
#****                                                           ****
#****                                                           ****
#****                                                           ****
#*******************************************************************
#*******************************************************************
DIR=${1%/*}/
NAME=${1#$DIR}
echo "DIR=$DIR"
echo "NAME=$NAME"
#====================================================================
#==  start writing parameter file  ber                             ==
#====================================================================
TMPFILE=${DIR}${NAME}.param
echo "write parameter file $TMPFILE"
echo "INPUT_DIR ${DIR}"               > ${TMPFILE}
echo "BASE_FILE_FORMAT PPM"           >> ${TMPFILE}
echo "INPUT_CONVERT $GIFTOPPM -verbose *" >> ${TMPFILE}
echo "INPUT"                          >> ${TMPFILE}
#=============================================================
#== list files here
#=============================================================
integer NIMAGE=0
DIGITS=
ZEROS=00000
for I in 1 2 3 4 5 6 7 8 9; do
  DIGITS=$DIGITS[0-9] 
  ZEROS=${ZEROS#0}
  for OLD in $DIR$NAME.${DIGITS}.gif ; do
    OLD=${OLD#$DIR}
    if [ -a $DIR$OLD ]; then
      let NIMAGE=NIMAGE+1
      echo "$OLD"  >>${TMPFILE}
    fi 
  done
done
echo "$NIMAGE gif files processd"
#=============================================================
#== all files listed
#=============================================================
echo "END_INPUT"                      >> ${TMPFILE}
echo "PATTERN IBBPBBPBBPBBPBB"        >> ${TMPFILE}
echo "FORCE_ENCODE_LAST_FRAME"        >> ${TMPFILE}
echo "GOP_SIZE 6"                     >> ${TMPFILE}
echo "SLICES_PER_FRAME 1"             >> ${TMPFILE}
echo "PIXEL HALF"                     >> ${TMPFILE}
echo "RANGE 10"                       >> ${TMPFILE}
echo "PSEARCH_ALG TWOLEVEL"           >> ${TMPFILE}
echo "BSEARCH_ALG CROSS2"             >> ${TMPFILE}
#*** This determines the quality/compression ratio
#*** lower numbers give better quality
echo "IQSCALE 8"                      >> ${TMPFILE}
echo "PQSCALE 10"                     >> ${TMPFILE}
echo "BQSCALE 11"                     >> ${TMPFILE}
echo "REFERENCE_FRAME DECODED"       >> ${TMPFILE}
echo "OUTPUT ${DIR}${NAME}.mpg"      >> ${TMPFILE}
#=============================================================
#== create mpg file 
#=============================================================
$MPEG_ENCODE ${TMPFILE}
#=============================================================
#== remove parameter file
#=============================================================

