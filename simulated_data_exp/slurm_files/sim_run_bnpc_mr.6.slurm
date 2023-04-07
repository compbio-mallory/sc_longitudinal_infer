#!/bin/bash
#SBATCH --job-name=bnpc6
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap001/rep2/m2/output6.out
#SBATCH --error=./sim_output/ap001/rep2/m2/error6.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap001/rep2/input_ap001_rep2.D.tsv -o ./sim_output/ap001/rep2/m2
