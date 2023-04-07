#!/bin/bash
#SBATCH --job-name=bnpc49
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap03/rep5/m5/output49.out
#SBATCH --error=./sim_output/ap03/rep5/m5/error49.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap03/rep5/input_ap03_rep5.D.tsv -o ./sim_output/ap03/rep5/m5
