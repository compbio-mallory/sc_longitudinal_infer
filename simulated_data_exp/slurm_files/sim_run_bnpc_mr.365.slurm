#!/bin/bash
#SBATCH --job-name=bnpc365
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep4/m1/output365.out
#SBATCH --error=./sim_output/theta3/rep4/m1/error365.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep4/input_theta3_rep4.D.tsv -o ./sim_output/theta3/rep4/m1
