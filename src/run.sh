# !/bin/bash
read -p "Enter the name of the system: " filename
declare -a BinCrossoverArray=("cpx" "onepx" "uniformx" "rrx" "mx")
declare -a IntCrossoverArray=("px1" "exchange")

eval "$(conda shell.bash hook)"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate general

for val in ${BinCrossoverArray[@]}; do
    for ix in ${IntCrossoverArray[@]}; do
        for ((i = 1; i<=31; i++)) do
            python3 main.py $i $filename $val $ix
            echo "Execution $i of $val and $ix finished"
        done
        echo "Evaluation of int crossover $ix Finished"
    done
    echo "Evaluation of $val Finished"
done

echo "Computation complete :)"
conda deactivate
