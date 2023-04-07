#!/bin/bash
#SBATCH --job-name=bnpc235
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mp3/rep3/m1/output235.out
#SBATCH --error=./sim_output/mp3/rep3/m1/error235.out
time python ../../BnpC/run_BnpC.py ./sim_input/mp3/rep3/input_mp3_rep3.D.tsv -o ./sim_output/mp3/rep3/m1
