#!/bin/bash
#SBATCH --job-name=bnpc158
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep2/m4/output158.out
#SBATCH --error=./sim_output/default/rep2/m4/error158.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep2/default_rep2.D.tsv -o ./sim_output/default/rep2/m4
