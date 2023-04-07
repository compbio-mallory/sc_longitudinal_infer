#!/bin/bash
#SBATCH --job-name=bnpc113
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep3/m4/output113.out
#SBATCH --error=./sim_output/bp3/rep3/m4/error113.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep3/input_bp3_rep3.D.tsv -o ./sim_output/bp3/rep3/m4
