# !/bin/bash
read -p "Enter the name of the system: " filename
read -p "Enter the binary crossover: " bincross
# read -p "Enter the theta value: " theta

declare -a IntCrossoverArray=("px1" "exchange")

eval "$(conda shell.bash hook)"
conda activate general

for ix in ${IntCrossoverArray[@]}; do
    for ((i = 1; i<= 31; i++)) do
	python3 main.py $i $filename $bincross $ix
	echo "Execution $i of $bincross and $ix finished"
    done
done
echo "Evaluation of $bincross complete :)"

conda deactivate
