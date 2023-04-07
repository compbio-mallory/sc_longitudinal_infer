#!/bin/bash
#SBATCH --job-name=bnpc151
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep1/m2/output151.out
#SBATCH --error=./sim_output/default/rep1/m2/error151.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep1/default_rep1.D.tsv -o ./sim_output/default/rep1/m2
