#!/bin/bash
#this script extracts the embedded source code from a cppaw binary
#call with: bash getsrc.sh $PAWBINARY $DESTINATION
#PAWBINARY is for example paw_fast.x
#DESTINATION is the location for the output of the source-tar-archive

PAW=$1
OUT=$2
start=`objdump -t $PAW | grep _binary_paw_srcblob_start | awk 'BEGIN { FS = " " } ; { print $1 }' | sed "s/^0*//g"`
end=`objdump -t $PAW | grep _binary_paw_srcblob_end | awk 'BEGIN { FS = " " } ; { print $1 }' | sed "s/^0*//g"`
size=`objdump -t $PAW | grep _binary_paw_srcblob_size | awk 'BEGIN { FS = " " } ; { print $1 }' | sed "s/^0*//g"`
startdata=`objdump -h $PAW | grep " .data " | awk 'BEGIN { FS = " " } ; { print $4 }' | sed "s/^0*//g"`
startdatareal=`objdump -h $PAW | grep " .data " | awk 'BEGIN { FS = " " } ; { print $6 }' | sed "s/^0*//g"`
startblobreal=`echo "$[0x$start]-$[0x$startdata]+$[0x$startdatareal]" | bc`
dd if=$PAW of=$OUT bs=1 count=$[0x$size] skip=$startblobreal
