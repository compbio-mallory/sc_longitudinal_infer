#!/bin/bash
#SBATCH --job-name=bnpc101
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp3/rep1/m2/output101.out
#SBATCH --error=./sim_output/bp3/rep1/m2/error101.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp3/rep1/input_bp3_rep1.D.tsv -o ./sim_output/bp3/rep1/m2
