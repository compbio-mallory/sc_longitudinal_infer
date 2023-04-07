#!/bin/bash
#SBATCH --job-name=bnpc168
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep4/m4/output168.out
#SBATCH --error=./sim_output/default/rep4/m4/error168.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep4/default_rep4.D.tsv -o ./sim_output/default/rep4/m4
