#!/bin/bash
#SBATCH --job-name=bnpc329
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep1/m5/output329.out
#SBATCH --error=./sim_output/t5/rep1/m5/error329.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep1/t5_rep1.D.tsv -o ./sim_output/t5/rep1/m5
