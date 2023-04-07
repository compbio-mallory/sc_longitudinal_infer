#!/bin/bash
#SBATCH --job-name=bnpc58
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap05/rep2/m4/output58.out
#SBATCH --error=./sim_output/ap05/rep2/m4/error58.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap05/rep2/input_ap05_rep2.D.tsv -o ./sim_output/ap05/rep2/m4
