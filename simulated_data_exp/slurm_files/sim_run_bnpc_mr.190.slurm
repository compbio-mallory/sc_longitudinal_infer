#!/bin/bash
#SBATCH --job-name=bnpc190
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep4/m1/output190.out
#SBATCH --error=./sim_output/mc12/rep4/m1/error190.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep4/input_mc12_rep4.D.tsv -o ./sim_output/mc12/rep4/m1
