#!/bin/bash
echo "comparing total energy."
#TOTAL_ENERGY_REF=-7.4079885 before commit 85f857388191c58f050acf5be749412fcfe54d12 (enforced that dtkin is symmetric) 
TOTAL_ENERGY_REF=-7.40813991 #tested with ifort12.1/MKL10.3 and gfortran 4.6.3/OpenBLAS/FFTW3

TOLERANCE=0.0001

TOTAL_ENERGY=`grep "TOTAL ENERGY" si2.prot | tail -n 1 | awk 'BEGIN { FS = " " } ; { print $4 }'`
CRIT=`echo "define abs(x) {if (x<0) {return -x}; return x;}scale=20;abs(($TOTAL_ENERGY_REF)-($TOTAL_ENERGY))<$TOLERANCE" | bc -l`
#echo $CRIT
if [ "$CRIT" = "1"  ];
then
  echo "TEST PASSED"
  exit 0
else
  echo "TEST FAILED"
  echo "SEE `pwd`."
  exit 1
fi
