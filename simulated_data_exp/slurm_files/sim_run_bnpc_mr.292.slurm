#!/bin/bash
#SBATCH --job-name=bnpc292
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep4/m3/output292.out
#SBATCH --error=./sim_output/sp5/rep4/m3/error292.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep4/input_sp5_rep4.D.tsv -o ./sim_output/sp5/rep4/m3
