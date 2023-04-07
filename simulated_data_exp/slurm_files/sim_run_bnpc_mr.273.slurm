#!/bin/bash
#SBATCH --job-name=bnpc273
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp05/rep5/m4/output273.out
#SBATCH --error=./sim_output/sp05/rep5/m4/error273.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp05/rep5/input_sp05_rep5.D.tsv -o ./sim_output/sp05/rep5/m4
