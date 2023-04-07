#!/bin/bash
#SBATCH --job-name=bnpc243
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mp3/rep4/m4/output243.out
#SBATCH --error=./sim_output/mp3/rep4/m4/error243.out
time python ../../BnpC/run_BnpC.py ./sim_input/mp3/rep4/input_mp3_rep4.D.tsv -o ./sim_output/mp3/rep4/m4
