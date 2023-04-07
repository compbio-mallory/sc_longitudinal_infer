#!/bin/bash
#SBATCH --job-name=bnpc146
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp4/rep5/m2/output146.out
#SBATCH --error=./sim_output/bp4/rep5/m2/error146.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp4/rep5/input_bp4_rep5.D.tsv -o ./sim_output/bp4/rep5/m2
