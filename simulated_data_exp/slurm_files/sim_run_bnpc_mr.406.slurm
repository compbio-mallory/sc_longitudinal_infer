#!/bin/bash
#SBATCH --job-name=bnpc406
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u0/rep2/m2/output406.out
#SBATCH --error=./sim_output/u0/rep2/m2/error406.out
time python ../../BnpC/run_BnpC.py ./sim_input/u0/rep2/input_u0_rep2.D.tsv -o ./sim_output/u0/rep2/m2
