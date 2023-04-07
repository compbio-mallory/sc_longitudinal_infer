#!/bin/bash
#SBATCH --job-name=bnpc159
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep2/m5/output159.out
#SBATCH --error=./sim_output/default/rep2/m5/error159.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep2/default_rep2.D.tsv -o ./sim_output/default/rep2/m5
