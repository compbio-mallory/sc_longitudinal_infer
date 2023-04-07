#!/bin/bash
#SBATCH --job-name=bnpc213
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep3/m4/output213.out
#SBATCH --error=./sim_output/mc2/rep3/m4/error213.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep3/input_mc2_rep3.D.tsv -o ./sim_output/mc2/rep3/m4
