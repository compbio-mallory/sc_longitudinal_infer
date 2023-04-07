#!/bin/bash
#SBATCH --job-name=bnpc249
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mp3/rep5/m5/output249.out
#SBATCH --error=./sim_output/mp3/rep5/m5/error249.out
time python ../../BnpC/run_BnpC.py ./sim_input/mp3/rep5/input_mp3_rep5.D.tsv -o ./sim_output/mp3/rep5/m5
