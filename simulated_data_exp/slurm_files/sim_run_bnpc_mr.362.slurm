#!/bin/bash
#SBATCH --job-name=bnpc362
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep3/m3/output362.out
#SBATCH --error=./sim_output/theta3/rep3/m3/error362.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep3/input_theta3_rep3.D.tsv -o ./sim_output/theta3/rep3/m3
