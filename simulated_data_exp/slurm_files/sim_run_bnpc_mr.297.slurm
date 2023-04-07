#!/bin/bash
#SBATCH --job-name=bnpc297
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep5/m3/output297.out
#SBATCH --error=./sim_output/sp5/rep5/m3/error297.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep5/input_sp5_rep5.D.tsv -o ./sim_output/sp5/rep5/m3
