#!/bin/bash
#SBATCH --job-name=bnpc29
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap03/rep1/m5/output29.out
#SBATCH --error=./sim_output/ap03/rep1/m5/error29.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap03/rep1/input_ap03_rep1.D.tsv -o ./sim_output/ap03/rep1/m5
