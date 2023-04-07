#!/bin/bash
#SBATCH --job-name=bnpc172
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep5/m3/output172.out
#SBATCH --error=./sim_output/default/rep5/m3/error172.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep5/default_rep5.D.tsv -o ./sim_output/default/rep5/m3
