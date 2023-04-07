#!/bin/bash
#SBATCH --job-name=bnpc301
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t4/rep1/m2/output301.out
#SBATCH --error=./sim_output/t4/rep1/m2/error301.out
time python ../../BnpC/run_BnpC.py ./sim_input/t4/rep1/t4_rep1.D.tsv -o ./sim_output/t4/rep1/m2
