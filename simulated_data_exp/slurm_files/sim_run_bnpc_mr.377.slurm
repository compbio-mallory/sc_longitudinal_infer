#!/bin/bash
#SBATCH --job-name=bnpc377
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep1/m3/output377.out
#SBATCH --error=./sim_output/theta5/rep1/m3/error377.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep1/input_theta5_rep1.D.tsv -o ./sim_output/theta5/rep1/m3
