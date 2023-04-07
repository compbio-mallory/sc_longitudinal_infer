#!/bin/bash
#SBATCH --job-name=bnpc209
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep2/m5/output209.out
#SBATCH --error=./sim_output/mc2/rep2/m5/error209.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep2/input_mc2_rep2.D.tsv -o ./sim_output/mc2/rep2/m5
