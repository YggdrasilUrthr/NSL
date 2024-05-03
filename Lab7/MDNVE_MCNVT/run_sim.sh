rm -rf output_*

phase="Solid"
verbose=false

while getopts slgv flag
do
    case "${flag}" in
        s) phase="Solid";;
        l) phase="Liquid";;
        g) phase="Gas";;
        v) verbose=true;;
    esac
done

echo "Executing simulation of $phase phase: "
cp ./NonEq/$phase/input.in .

if $verbose ; then
    ./NVE_NVT.out
else
    echo "Running equilbration step..."
    ./NVE_NVT.out >/dev/null
fi

rm -rf ./output_*
cp ./Eq/$phase/input.in .

if $verbose ; then
    ./NVE_NVT.out
else
    echo "Running simulation step..."
    ./NVE_NVT.out >/dev/null
fi

mv ./output_* ./Eq/$phase/
echo "Done"