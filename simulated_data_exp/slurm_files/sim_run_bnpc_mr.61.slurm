#!/bin/bash
#SBATCH --job-name=bnpc61
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap05/rep3/m2/output61.out
#SBATCH --error=./sim_output/ap05/rep3/m2/error61.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap05/rep3/input_ap05_rep3.D.tsv -o ./sim_output/ap05/rep3/m2
