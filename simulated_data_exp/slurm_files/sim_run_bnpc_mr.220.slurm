#!/bin/bash
#SBATCH --job-name=bnpc220
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep5/m1/output220.out
#SBATCH --error=./sim_output/mc2/rep5/m1/error220.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep5/input_mc2_rep5.D.tsv -o ./sim_output/mc2/rep5/m1
