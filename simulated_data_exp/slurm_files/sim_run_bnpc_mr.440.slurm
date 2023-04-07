#!/bin/bash
#SBATCH --job-name=bnpc440
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u1/rep4/m1/output440.out
#SBATCH --error=./sim_output/u1/rep4/m1/error440.out
time python ../../BnpC/run_BnpC.py ./sim_input/u1/rep4/input_u1_rep4.D.tsv -o ./sim_output/u1/rep4/m1
