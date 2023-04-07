#!/bin/bash
#SBATCH --job-name=bnpc450
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep1/m1/output450.out
#SBATCH --error=./sim_output/u3/rep1/m1/error450.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep1/input_u3_rep1.D.tsv -o ./sim_output/u3/rep1/m1
