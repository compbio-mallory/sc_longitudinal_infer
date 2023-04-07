#!/bin/bash
#SBATCH --job-name=bnpc424
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u0/rep5/m5/output424.out
#SBATCH --error=./sim_output/u0/rep5/m5/error424.out
time python ../../BnpC/run_BnpC.py ./sim_input/u0/rep5/input_u0_rep5.D.tsv -o ./sim_output/u0/rep5/m5
