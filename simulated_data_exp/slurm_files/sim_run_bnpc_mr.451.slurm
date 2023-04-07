#!/bin/bash
#SBATCH --job-name=bnpc451
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep1/m2/output451.out
#SBATCH --error=./sim_output/u3/rep1/m2/error451.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep1/input_u3_rep1.D.tsv -o ./sim_output/u3/rep1/m2
