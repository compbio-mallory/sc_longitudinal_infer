#!/bin/bash
#SBATCH --job-name=bnpc459
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep2/m5/output459.out
#SBATCH --error=./sim_output/u3/rep2/m5/error459.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep2/input_u3_rep2.D.tsv -o ./sim_output/u3/rep2/m5
