#!/bin/bash
#SBATCH --job-name=bnpc170
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/default/rep5/m1/output170.out
#SBATCH --error=./sim_output/default/rep5/m1/error170.out
time python ../../BnpC/run_BnpC.py ./sim_input/default/rep5/default_rep5.D.tsv -o ./sim_output/default/rep5/m1
