#!/bin/bash
#SBATCH --job-name=bnpc344
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep4/m5/output344.out
#SBATCH --error=./sim_output/t5/rep4/m5/error344.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep4/t5_rep4.D.tsv -o ./sim_output/t5/rep4/m5
