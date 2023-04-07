#!/bin/bash
#SBATCH --job-name=bnpc426
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u1/rep1/m2/output426.out
#SBATCH --error=./sim_output/u1/rep1/m2/error426.out
time python ../../BnpC/run_BnpC.py ./sim_input/u1/rep1/input_u1_rep1.D.tsv -o ./sim_output/u1/rep1/m2
