#!/bin/bash

# SPFUNCTION_BESSEL "sqrt(Pi/(2*x))*BesselJ(l+1/2,x)"
# SPFUNCTION_BESSEL0 "x^l/product(2*i-1,i=1..l+1)"
# SPFUNCTION_NEUMANN "sqrt(Pi/(2*x))*BesselY(l+1/2,x)"
# SPFUNCTION_NEUMANN0 "-product(2*i-1,i=1..l)/x^(l+1)"
# SPFUNCTION_MODBESSEL "sqrt(Pi/(2*x))*BesselI(l+1/2,x)"
# SPFUNCTION_MODNEUMANN "sqrt(Pi/(2*x))*BesselI(-l-1/2,x)"
# SPFUNCTION_MODHANKEL "Pi/2*sqrt(Pi/(2*x))*(-1)^(l+1)*(BesselI(l+1/2,x)-BesselI(-l-1/2,x))"

rm refdata.dat
MAPLE="/opt/maple17/bin/maple"
findex="0"
for function in "sqrt(Pi/(2*x))*BesselJ(l+1/2,x)" "x^l/product(2*i-1,i=1..l+1)" "sqrt(Pi/(2*x))*BesselY(l+1/2,x)" "-product(2*i-1,i=1..l)/x^(l+1)" "sqrt(Pi/(2*x))*BesselI(l+1/2,x)" "sqrt(Pi/(2*x))*BesselI(-l-1/2,x)" "Pi/2*sqrt(Pi/(2*x))*(-1)^(l+1)*(BesselI(l+1/2,x)-BesselI(-l-1/2,x))";
do
  findex=`echo "($findex+1)" | bc -l`
  for order in `seq 0 20`;
  do
    echo "$function $order"
    xindex="0"
    for logx in `seq -10 6`;
    do
      xindex=`echo "($xindex+1)" | bc -l`
      Y=`echo "x:=10^($logx);l:=$order;y:=evalf( $function ,500);printf(\"z%.30g\",y);" | $MAPLE -q | grep "z"| sed 's/[^-.0-9eE]//g' | tr -d '\n'  | sed "s/e/D/g"`
      DYDX=`echo "X:=10^($logx);l:=$order;y:=evalf(limit(diff( $function ,x),x=X),500);printf(\"z%.30g\",y);" | $MAPLE -q | grep "z"| sed 's/[^-.0-9eE]//g' | tr -d '\n' | sed "s/e/D/g"`
      hasD=`echo "$Y" | grep "D" | wc -l`
      if [ "$hasD" = "0" ];
      then
        y="${Y}D0"
        Y="$y"
      fi
      hasD=`echo "$DYDX" | grep "D" | wc -l`
      if [ "$hasD" = "0" ];
      then
        dydx="${DYDX}D0"
        DYDX="$dydx"
      fi
      Y2=`echo "define abs(x) { if ( x<0 ) return -x; return x };scale=2000;abs($Y))>10^300" | sed "s/D/\*10\^(/g" | bc -l`
      DYDX2=`echo "define abs(x) { if ( x<0 ) return -x; return x };scale=2000;abs($DYDX))>10^300" | sed "s/D/\*10\^(/g" | bc -l`
      Y3=`echo "define abs(x) { if ( x<0 ) return -x; return x };scale=2000;abs($Y))<10^(-300)" | sed "s/D/\*10\^(/g" | bc -l`
      DYDX3=`echo "define abs(x) { if ( x<0 ) return -x; return x };scale=2000;abs($DYDX))<10^(-300)" | sed "s/D/\*10\^(/g" | bc -l`
      err=`echo "$Y2+$DYDX2+$Y3+$DYDX3" | bc -l`
      if [ "$err" = "0" ];
      then  
          echo "REFDATA_Y($findex,$order,$xindex)=$Y" >> refdata.dat
          echo "REFDATA_DYDX($findex,$order,$xindex)=$DYDX" >> refdata.dat
      else
          echo "REFDATA_Y($findex,$order,$xindex)=0.0D0! $Y" >> refdata.dat
          echo "REFDATA_DYDX($findex,$order,$xindex)=0.0D0! $DYDX" >> refdata.dat
          echo "EXCLUDE($findex,$order,$xindex)=.true." >> refdata.dat

      fi
    done
  done
done


#y=`echo "x:=$x;evalf(sqrt(Pi/(2*x))*BesselJ($order+1/2,x),200);" | /opt/maple16/bin/maple -q | grep -v "x"| sed 's/[^.0-9]//g'`
#echo "$x\t$order\t$y"dat
