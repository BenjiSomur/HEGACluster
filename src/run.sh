# !/bin/bash
read -p "Enter the name of the system: " filename
read -p "Enter the theta value: " theta
# read -p "Enter the crossover method: " copx
# read -p "Enter the int crossover method: " intx
declare -a BinCrossoverArray=("cpx" "onepx")
declare -a IntCrossoverArray=("px1" "exchange")
eval "$(conda shell.bash hook)"
conda activate general

for val in ${BinCrossoverArray[@]}; do
    for ix in ${IntCrossoverArray[@]}; do
        for ((i = 1; i<=31; i++)) do
            python3 main.py $i $filename $theta $val $ix
            echo "Execution $i of $val and $ix finished"
        done
        echo "Evaluation of int crossover $ix Finished"
    done
    echo "Evaluation of bin crossover $val Finished"
done

echo "Computation complete :)"
conda deactivate
