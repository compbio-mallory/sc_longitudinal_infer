#!/bin/bash
#SBATCH --job-name=bnpc380
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep2/m1/output380.out
#SBATCH --error=./sim_output/theta5/rep2/m1/error380.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep2/input_theta5_rep2.D.tsv -o ./sim_output/theta5/rep2/m1
