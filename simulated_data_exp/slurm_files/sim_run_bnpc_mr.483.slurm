#!/bin/bash
#SBATCH --job-name=bnpc483
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/v1v2/rep2/m4/output483.out
#SBATCH --error=./sim_output/v1v2/rep2/m4/error483.out
time python ../../BnpC/run_BnpC.py ./sim_input/v1v2/rep2/input_v1v2_rep2.D.tsv -o ./sim_output/v1v2/rep2/m4
