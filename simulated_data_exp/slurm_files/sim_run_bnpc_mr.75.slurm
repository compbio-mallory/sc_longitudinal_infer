#!/bin/bash
#SBATCH --job-name=bnpc75
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp1/rep1/m1/output75.out
#SBATCH --error=./sim_output/bp1/rep1/m1/error75.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp1/rep1/input_bp1_rep1.D.tsv -o ./sim_output/bp1/rep1/m1
