#!/bin/bash
#SBATCH --job-name=bnpc216
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep4/m2/output216.out
#SBATCH --error=./sim_output/mc2/rep4/m2/error216.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep4/input_mc2_rep4.D.tsv -o ./sim_output/mc2/rep4/m2
