#!/bin/bash

TStart=0.5
TStop=2.0
NStep=20

h=${1:-0}
m=${2:-1}
g=""

if [ $m == "0" ]; then
    g="_g"
fi

cat /dev/null > output.chi.$h$g
cat /dev/null > output.ene.$h$g
cat /dev/null > output.hec.$h$g
cat /dev/null > output.mag.$h$g

sed -i '4s/.*/'$h'/' input.dat
sed -i '5s/.*/'$m'/' input.dat

for i in {0..20}
do
#    T=$(echo "scale=2; $i*($TStop-$TStart)/$NStep+$TStart" | bc -l)
    T=$(echo "scale=2; ($NStep - $i)*($TStop-$TStart)/$NStep+$TStart" | bc -l)
    sed -i '1s/.*/'$T'/' input.dat
    ./Monte_Carlo_ISING_1D.out
done
