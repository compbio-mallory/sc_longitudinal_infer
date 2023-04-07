#!/bin/bash
#SBATCH --job-name=bnpc328
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep1/m4/output328.out
#SBATCH --error=./sim_output/t5/rep1/m4/error328.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep1/t5_rep1.D.tsv -o ./sim_output/t5/rep1/m4
