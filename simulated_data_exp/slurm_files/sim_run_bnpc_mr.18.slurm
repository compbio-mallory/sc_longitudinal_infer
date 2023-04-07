#!/bin/bash
#SBATCH --job-name=bnpc18
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap001/rep4/m4/output18.out
#SBATCH --error=./sim_output/ap001/rep4/m4/error18.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap001/rep4/input_ap001_rep4.D.tsv -o ./sim_output/ap001/rep4/m4
