#!/bin/sh
# 
# copies all files ending with string $1 into files beginning with
# the same beginning and ending in $2
#
for FILE in *$1 ;do cp -p "$FILE" "${FILE%$1}$2"; done
exit
