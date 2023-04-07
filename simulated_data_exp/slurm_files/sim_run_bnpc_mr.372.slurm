#!/bin/bash
#SBATCH --job-name=bnpc372
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep5/m3/output372.out
#SBATCH --error=./sim_output/theta3/rep5/m3/error372.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep5/input_theta3_rep5.D.tsv -o ./sim_output/theta3/rep5/m3
