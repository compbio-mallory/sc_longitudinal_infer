#!/bin/bash
#SBATCH --job-name=bnpc87
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp1/rep3/m3/output87.out
#SBATCH --error=./sim_output/bp1/rep3/m3/error87.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp1/rep3/input_bp1_rep3.D.tsv -o ./sim_output/bp1/rep3/m3
