#!/bin/bash
echo "comparing total energy."
TOTAL_ENERGY_REF=-7.4079885
TOLERANCE=0.0001

TOTAL_ENERGY=`grep "TOTAL ENERGY" si2.prot | tail -n 1 | awk 'BEGIN { FS = " " } ; { print $4 }'`
CRIT=`echo "define abs(x) {if (x<0) {return -x}; return x;}scale=20;abs(($TOTAL_ENERGY_REF)-($TOTAL_ENERGY))<$TOLERANCE" | bc -l`
echo $CRIT
if [ "$CRIT" = "1"  ];
then
  echo "TEST PASSED"
  exit 0
else
  echo "TEST FAILED"
  echo "SEE `pwd`."
  exit 1
fi
