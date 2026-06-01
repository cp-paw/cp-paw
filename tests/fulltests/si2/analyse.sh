#!/bin/bash
echo "comparing total energy."
TOTAL_ENERGY_REF=-7.3657485
# PBE TYPE=10 after the Dec. 10, 2023 PBE TFAC fix. The legacy TYPE=-10 path
# reproduces the previous reference, -7.4076060.

TOLERANCE=0.0001

TOTAL_ENERGY=`grep "TOTAL ENERGY" si2.prot | tail -n 1 | awk 'BEGIN { FS = " " } ; { print $4 }'`
CRIT=`awk -v ref="$TOTAL_ENERGY_REF" -v val="$TOTAL_ENERGY" -v tol="$TOLERANCE" \
  'BEGIN { diff=ref-val; if (diff < 0) diff=-diff; print (diff < tol) ? 1 : 0 }'`
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
