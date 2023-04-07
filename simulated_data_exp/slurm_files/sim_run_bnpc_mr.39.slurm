#!/bin/bash
#SBATCH --job-name=bnpc39
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap03/rep3/m5/output39.out
#SBATCH --error=./sim_output/ap03/rep3/m5/error39.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap03/rep3/input_ap03_rep3.D.tsv -o ./sim_output/ap03/rep3/m5
