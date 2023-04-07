#!/bin/bash
#SBATCH --job-name=bnpc447
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u1/rep5/m3/output447.out
#SBATCH --error=./sim_output/u1/rep5/m3/error447.out
time python ../../BnpC/run_BnpC.py ./sim_input/u1/rep5/input_u1_rep5.D.tsv -o ./sim_output/u1/rep5/m3
