#!/bin/bash
#SBATCH --job-name=bnpc491
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/v1v2/rep4/m2/output491.out
#SBATCH --error=./sim_output/v1v2/rep4/m2/error491.out
time python ../../BnpC/run_BnpC.py ./sim_input/v1v2/rep4/input_v1v2_rep4.D.tsv -o ./sim_output/v1v2/rep4/m2
