#!/bin/bash
#SBATCH --job-name=bnpc286
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep3/m2/output286.out
#SBATCH --error=./sim_output/sp5/rep3/m2/error286.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep3/input_sp5_rep3.D.tsv -o ./sim_output/sp5/rep3/m2
