#!/bin/bash
#SBATCH --job-name=bnpc198
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep5/m4/output198.out
#SBATCH --error=./sim_output/mc12/rep5/m4/error198.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep5/input_mc12_rep5.D.tsv -o ./sim_output/mc12/rep5/m4
