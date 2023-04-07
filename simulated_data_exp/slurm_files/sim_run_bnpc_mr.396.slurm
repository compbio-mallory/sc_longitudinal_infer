#!/bin/bash
#SBATCH --job-name=bnpc396
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep5/m2/output396.out
#SBATCH --error=./sim_output/theta5/rep5/m2/error396.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep5/input_theta5_rep5.D.tsv -o ./sim_output/theta5/rep5/m2
