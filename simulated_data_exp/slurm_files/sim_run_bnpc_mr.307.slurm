#!/bin/bash
#SBATCH --job-name=bnpc307
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t4/rep2/m3/output307.out
#SBATCH --error=./sim_output/t4/rep2/m3/error307.out
time python ../../BnpC/run_BnpC.py ./sim_input/t4/rep2/t4_rep2.D.tsv -o ./sim_output/t4/rep2/m3
