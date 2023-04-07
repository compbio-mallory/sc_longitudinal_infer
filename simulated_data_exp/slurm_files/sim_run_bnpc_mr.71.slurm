#!/bin/bash
#SBATCH --job-name=bnpc71
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap05/rep5/m2/output71.out
#SBATCH --error=./sim_output/ap05/rep5/m2/error71.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap05/rep5/input_ap05_rep5.D.tsv -o ./sim_output/ap05/rep5/m2
