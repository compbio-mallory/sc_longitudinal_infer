#!/bin/bash
#SBATCH --job-name=bnpc392
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep4/m3/output392.out
#SBATCH --error=./sim_output/theta5/rep4/m3/error392.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep4/input_theta5_rep4.D.tsv -o ./sim_output/theta5/rep4/m3
