#!/bin/bash
#SBATCH --job-name=bnpc333
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep2/m4/output333.out
#SBATCH --error=./sim_output/t5/rep2/m4/error333.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep2/t5_rep2.D.tsv -o ./sim_output/t5/rep2/m4
