#!/bin/bash
#SBATCH --job-name=bnpc351
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep1/m2/output351.out
#SBATCH --error=./sim_output/theta3/rep1/m2/error351.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep1/input_theta3_rep1.D.tsv -o ./sim_output/theta3/rep1/m2
