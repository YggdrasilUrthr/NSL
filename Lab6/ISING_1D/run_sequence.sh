#!/bin/bash

TStart=0.5
TStop=2.0
NStep=20

cat /dev/null > output.chi.0
cat /dev/null > output.ene.0
cat /dev/null > output.mag.0

for i in {0..20}
do
    T=$(echo "scale=2; $i*($TStop-$TStart)/$NStep+$TStart" | bc -l)
    sed -i '1s/.*/'$T'/' input.dat
    ./Monte_Carlo_ISING_1D.out
done
