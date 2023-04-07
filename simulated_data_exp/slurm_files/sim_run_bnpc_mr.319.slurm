#!/bin/bash
#SBATCH --job-name=bnpc319
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t4/rep4/m5/output319.out
#SBATCH --error=./sim_output/t4/rep4/m5/error319.out
time python ../../BnpC/run_BnpC.py ./sim_input/t4/rep4/t4_rep4.D.tsv -o ./sim_output/t4/rep4/m5
