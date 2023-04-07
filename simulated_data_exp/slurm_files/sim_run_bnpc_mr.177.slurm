#!/bin/bash
#SBATCH --job-name=bnpc177
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep1/m3/output177.out
#SBATCH --error=./sim_output/mc12/rep1/m3/error177.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep1/input_mc12_rep1.D.tsv -o ./sim_output/mc12/rep1/m3
