#!/bin/bash
#SBATCH --job-name=bnpc254
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp05/rep1/m5/output254.out
#SBATCH --error=./sim_output/sp05/rep1/m5/error254.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp05/rep1/input_sp05_rep1.D.tsv -o ./sim_output/sp05/rep1/m5
