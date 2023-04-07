#!/bin/bash
#SBATCH --job-name=bnpc300
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t4/rep1/m1/output300.out
#SBATCH --error=./sim_output/t4/rep1/m1/error300.out
time python ../../BnpC/run_BnpC.py ./sim_input/t4/rep1/t4_rep1.D.tsv -o ./sim_output/t4/rep1/m1
