#!/bin/bash
#SBATCH --job-name=bnpc267
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp05/rep4/m3/output267.out
#SBATCH --error=./sim_output/sp05/rep4/m3/error267.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp05/rep4/input_sp05_rep4.D.tsv -o ./sim_output/sp05/rep4/m3
