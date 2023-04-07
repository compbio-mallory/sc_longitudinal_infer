#!/bin/bash
#SBATCH --job-name=bnpc107
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep2/m3/output107.out
#SBATCH --error=./sim_output/bp3/rep2/m3/error107.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep2/input_bp3_rep2.D.tsv -o ./sim_output/bp3/rep2/m3
