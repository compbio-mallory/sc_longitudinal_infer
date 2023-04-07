#!/bin/bash
#SBATCH --job-name=bnpc256
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp05/rep2/m2/output256.out
#SBATCH --error=./sim_output/sp05/rep2/m2/error256.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp05/rep2/input_sp05_rep2.D.tsv -o ./sim_output/sp05/rep2/m2
