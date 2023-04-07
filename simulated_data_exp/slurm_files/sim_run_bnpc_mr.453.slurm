#!/bin/bash
#SBATCH --job-name=bnpc453
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep1/m4/output453.out
#SBATCH --error=./sim_output/u3/rep1/m4/error453.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep1/input_u3_rep1.D.tsv -o ./sim_output/u3/rep1/m4
