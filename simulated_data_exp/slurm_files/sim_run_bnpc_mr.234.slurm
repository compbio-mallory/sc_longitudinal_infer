#!/bin/bash
#SBATCH --job-name=bnpc234
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mp3/rep2/m5/output234.out
#SBATCH --error=./sim_output/mp3/rep2/m5/error234.out
time python ../../BnpC/run_BnpC.py ./sim_input/mp3/rep2/input_mp3_rep2.D.tsv -o ./sim_output/mp3/rep2/m5
