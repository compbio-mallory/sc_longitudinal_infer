#!/bin/bash
#SBATCH --job-name=bnpc121
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep5/m2/output121.out
#SBATCH --error=./sim_output/bp3/rep5/m2/error121.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep5/input_bp3_rep5.D.tsv -o ./sim_output/bp3/rep5/m2
