rm -rf output_*

simpth="./ArSim/"
phase="Solid"
smet="MD"
verbose=false

while getopts slgmv flag
do
    case "${flag}" in
        s) phase="Solid";;
        l) phase="Liquid";;
        g) phase="Gas";;
        m) smet="MC";;          # Select montecarlo sampling
        v) verbose=true;;
    esac
done

echo "Executing simulation of Argon $phase phase, using $smet sampling: "
cp $simpth$smet/NonEq/$phase/input.in .

if $verbose ; then
    ./NVE_NVT.out
else
    echo "Running equilbration step..."
    ./NVE_NVT.out >/dev/null
fi

rm -rf ./output_*
cp $simpth$smet/Eq/$phase/input.in .

if $verbose ; then
    ./NVE_NVT.out
else
    echo "Running simulation step..."
    ./NVE_NVT.out >/dev/null
fi

mv ./output_* $simpth$smet/Eq/$phase/
echo "Done"