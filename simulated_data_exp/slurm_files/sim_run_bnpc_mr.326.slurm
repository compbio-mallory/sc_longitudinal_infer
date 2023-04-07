#!/bin/bash
#SBATCH --job-name=bnpc326
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep1/m2/output326.out
#SBATCH --error=./sim_output/t5/rep1/m2/error326.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep1/t5_rep1.D.tsv -o ./sim_output/t5/rep1/m2
