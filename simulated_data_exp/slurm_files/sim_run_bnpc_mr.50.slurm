#!/bin/bash
#SBATCH --job-name=bnpc50
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap05/rep1/m1/output50.out
#SBATCH --error=./sim_output/ap05/rep1/m1/error50.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap05/rep1/input_ap05_rep1.D.tsv -o ./sim_output/ap05/rep1/m1
