#!/bin/bash
#SBATCH --job-name=bnpc171
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep5/m2/output171.out
#SBATCH --error=./sim_output/default/rep5/m2/error171.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep5/default_rep5.D.tsv -o ./sim_output/default/rep5/m2
