#!/bin/bash
#SBATCH --job-name=bnpc439
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u1/rep3/m5/output439.out
#SBATCH --error=./sim_output/u1/rep3/m5/error439.out
time python ../../BnpC/run_BnpC.py ./sim_input/u1/rep3/input_u1_rep3.D.tsv -o ./sim_output/u1/rep3/m5
