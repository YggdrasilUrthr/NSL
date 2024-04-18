rm -rf output_*

phase="Solid"

while getopts slg flag
do
    case "${flag}" in
        s) phase="Solid";;
        l) phase="Liquid";;
        g) phase="Gas";;
    esac
done

cp ./NonEq/$phase/input.in .
./NVE_NVT.out
cp ./output_temp.dat ./NonEq/$phase/