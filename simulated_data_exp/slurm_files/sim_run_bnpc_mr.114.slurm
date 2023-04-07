#!/bin/bash
#SBATCH --job-name=bnpc114
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep3/m5/output114.out
#SBATCH --error=./sim_output/bp3/rep3/m5/error114.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep3/input_bp3_rep3.D.tsv -o ./sim_output/bp3/rep3/m5
