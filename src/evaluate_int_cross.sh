# !/bin/bash
read -p "Enter the name of the system: " filename
read -p "Enter the int crossover: " intcross
read -p "Enter the theta value: " theta

declare -a BinCrossoverArray=("cpx onepx")

eval "$(conda shell.bash hook)"
conda activate general

for val in ${BinCrossoverArray[@]}; do
    for ((i=1; i<=31; i++)) do
        python3 maint.py $i $filename $theta $val $intcross
        echo "Execution $i of $val and $intcross finished"
    done
done

echo "Evaluation of $intcross complete :)"
conda deactivate

