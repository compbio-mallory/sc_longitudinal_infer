#!/bin/bash
#SBATCH --job-name=bnpc489
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/v1v2/rep3/m5/output489.out
#SBATCH --error=./sim_output/v1v2/rep3/m5/error489.out
time python ../../BnpC/run_BnpC.py ./sim_input/v1v2/rep3/input_v1v2_rep3.D.tsv -o ./sim_output/v1v2/rep3/m5
