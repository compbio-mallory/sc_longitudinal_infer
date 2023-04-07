#!/bin/bash
#SBATCH --job-name=bnpc33
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap03/rep2/m4/output33.out
#SBATCH --error=./sim_output/ap03/rep2/m4/error33.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap03/rep2/input_ap03_rep2.D.tsv -o ./sim_output/ap03/rep2/m4
