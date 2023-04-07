#!/bin/bash
#SBATCH --job-name=bnpc225
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mp3/rep1/m1/output225.out
#SBATCH --error=./sim_output/mp3/rep1/m1/error225.out
time python ../../BnpC/run_BnpC.py ./sim_input/mp3/rep1/input_mp3_rep1.D.tsv -o ./sim_output/mp3/rep1/m1
