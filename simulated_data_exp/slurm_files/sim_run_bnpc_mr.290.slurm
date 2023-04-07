#!/bin/bash
#SBATCH --job-name=bnpc290
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep4/m1/output290.out
#SBATCH --error=./sim_output/sp5/rep4/m1/error290.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep4/input_sp5_rep4.D.tsv -o ./sim_output/sp5/rep4/m1
