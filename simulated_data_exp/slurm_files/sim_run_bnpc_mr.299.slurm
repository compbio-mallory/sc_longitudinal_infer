#!/bin/bash
#SBATCH --job-name=bnpc299
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep5/m5/output299.out
#SBATCH --error=./sim_output/sp5/rep5/m5/error299.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep5/input_sp5_rep5.D.tsv -o ./sim_output/sp5/rep5/m5
