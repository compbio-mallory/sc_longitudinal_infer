#!/bin/bash
#SBATCH --job-name=bnpc52
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap05/rep1/m3/output52.out
#SBATCH --error=./sim_output/ap05/rep1/m3/error52.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap05/rep1/input_ap05_rep1.D.tsv -o ./sim_output/ap05/rep1/m3
