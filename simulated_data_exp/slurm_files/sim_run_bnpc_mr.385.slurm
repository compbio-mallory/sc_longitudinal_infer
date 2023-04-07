#!/bin/bash
#SBATCH --job-name=bnpc385
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep3/m1/output385.out
#SBATCH --error=./sim_output/theta5/rep3/m1/error385.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep3/input_theta5_rep3.D.tsv -o ./sim_output/theta5/rep3/m1
