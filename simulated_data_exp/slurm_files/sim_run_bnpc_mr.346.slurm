#!/bin/bash
#SBATCH --job-name=bnpc346
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep5/m2/output346.out
#SBATCH --error=./sim_output/t5/rep5/m2/error346.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep5/t5_rep5.D.tsv -o ./sim_output/t5/rep5/m2
