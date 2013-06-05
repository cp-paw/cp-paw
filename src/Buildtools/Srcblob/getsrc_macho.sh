#!/bin/bash
#this script extracts the embedded source code from a cppaw binary which is in macho-format
#if you have a binary that has been created on linux/BSD please use getsrc_elf.sh
#call with: bash getsrc_elf.sh $PAWBINARY $DESTINATION
#PAWBINARY is for example paw_fast.x
#DESTINATION is the location for the output of the source-tar-archive

PAW=$1
OUT=$2
otool $PAW -s binary pawsrcblob_bin > pawsrc_tmp
tail +3 pawsrc_tmp | awk 'BEGIN { FS = " " } ; { print $2 $3 $4 $5 $6 $7 $8 $9 $10 $11 $12 $13 $14 $15 $16 $17 }' > pawsrc_tmp2
xxd -r -p pawsrc_tmp2 $OUT
rm pawsrc_tmp
rm pawsrc_tmp2
