#!/bin/bash

# Runs a simulation sequence at decreasing temperatures, dumping equilibration data if needed

TStart=0.5
TStop=2.0
NStep=20

CSV="./CSV/"
SMet="Metr/"
Eq="Eq/"

h=${1:-0}       # h parameter, default 0
m=${2:-1}       # sampling method, default Metropolis 
e=${@: 3}       # dump equilibration data at specified temperatures

g=""

if [ $m == "0" ]; then
    g="_g"
    SMet="Gibbs/"
fi

rm -rf output.*

for i in {0..20}
do
    # Use 4 decimal places to be on the safe side (2 gave me problems for certain values of NStep)
    T=$(echo "scale=4; ($NStep - $i)*($TStop-$TStart)/$NStep+$TStart" | bc -l)
    echo $T

    # Equilibration
    cp ./Input/input_eq.dat ./input.dat
    sed -i '1s/.*/'$T'/' input.dat
    sed -i '4s/.*/'$h'/' input.dat
    sed -i '5s/.*/'$m'/' input.dat
    ./Monte_Carlo_ISING_1D.out
    if [[ ${e[@]} =~ $T ]]; then
        for f in eq_output.mag.*; do mv $f $CSV$Eq$SMet$f"_"$T; done
    else
        # These can be safely removed, only config.final is needed
        rm -rf eq_output.*
    fi
    
    # Simulation
    cp ./Input/input_run.dat ./input.dat
    sed -i '1s/.*/'$T'/' input.dat
    sed -i '4s/.*/'$h'/' input.dat
    sed -i '5s/.*/'$m'/' input.dat
    ./Monte_Carlo_ISING_1D.out
done

mv output.* $CSV$SMet
