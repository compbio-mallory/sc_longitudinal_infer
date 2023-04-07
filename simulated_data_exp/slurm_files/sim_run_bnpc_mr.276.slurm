#!/bin/bash
#SBATCH --job-name=bnpc276
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep1/m2/output276.out
#SBATCH --error=./sim_output/sp5/rep1/m2/error276.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep1/input_sp5_rep1.D.tsv -o ./sim_output/sp5/rep1/m2
