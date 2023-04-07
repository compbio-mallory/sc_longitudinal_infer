#!/bin/bash
#SBATCH --job-name=bnpc116
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep4/m2/output116.out
#SBATCH --error=./sim_output/bp3/rep4/m2/error116.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep4/input_bp3_rep4.D.tsv -o ./sim_output/bp3/rep4/m2
