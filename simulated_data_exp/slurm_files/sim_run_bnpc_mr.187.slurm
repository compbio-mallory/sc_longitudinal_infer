#!/bin/bash
#SBATCH --job-name=bnpc187
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep3/m3/output187.out
#SBATCH --error=./sim_output/mc12/rep3/m3/error187.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep3/input_mc12_rep3.D.tsv -o ./sim_output/mc12/rep3/m3
