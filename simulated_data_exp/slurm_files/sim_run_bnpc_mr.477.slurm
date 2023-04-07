#!/bin/bash
#SBATCH --job-name=bnpc477
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/v1v2/rep1/m3/output477.out
#SBATCH --error=./sim_output/v1v2/rep1/m3/error477.out
time python ../../BnpC/run_BnpC.py ./sim_input/v1v2/rep1/input_v1v2_rep1.D.tsv -o ./sim_output/v1v2/rep1/m3
