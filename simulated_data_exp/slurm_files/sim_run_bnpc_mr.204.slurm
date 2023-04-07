#!/bin/bash
#SBATCH --job-name=bnpc204
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep1/m5/output204.out
#SBATCH --error=./sim_output/mc2/rep1/m5/error204.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep1/input_mc2_rep1.D.tsv -o ./sim_output/mc2/rep1/m5
