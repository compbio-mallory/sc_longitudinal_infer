#!/bin/bash
#SBATCH --job-name=bnpc97
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp1/rep5/m3/output97.out
#SBATCH --error=./sim_output/bp1/rep5/m3/error97.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp1/rep5/input_bp1_rep5.D.tsv -o ./sim_output/bp1/rep5/m3
