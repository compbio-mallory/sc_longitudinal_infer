#!/bin/bash
#SBATCH --job-name=bnpc161
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep3/m2/output161.out
#SBATCH --error=./sim_output/default/rep3/m2/error161.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep3/default_rep3.D.tsv -o ./sim_output/default/rep3/m2
