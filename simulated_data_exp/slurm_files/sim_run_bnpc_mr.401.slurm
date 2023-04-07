#!/bin/bash
#SBATCH --job-name=bnpc401
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u0/rep1/m2/output401.out
#SBATCH --error=./sim_output/u0/rep1/m2/error401.out
time python ../../BnpC/run_BnpC.py ./sim_input/u0/rep1/input_u0_rep1.D.tsv -o ./sim_output/u0/rep1/m2
