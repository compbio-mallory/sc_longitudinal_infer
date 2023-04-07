#!/bin/bash
#SBATCH --job-name=bnpc334
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep2/m5/output334.out
#SBATCH --error=./sim_output/t5/rep2/m5/error334.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep2/t5_rep2.D.tsv -o ./sim_output/t5/rep2/m5
