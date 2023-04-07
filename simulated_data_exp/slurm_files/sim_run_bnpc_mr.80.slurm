#!/bin/bash
#SBATCH --job-name=bnpc80
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp1/rep2/m1/output80.out
#SBATCH --error=./sim_output/bp1/rep2/m1/error80.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp1/rep2/input_bp1_rep2.D.tsv -o ./sim_output/bp1/rep2/m1
