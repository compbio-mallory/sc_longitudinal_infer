#!/bin/bash
#SBATCH --job-name=bnpc22
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap001/rep5/m3/output22.out
#SBATCH --error=./sim_output/ap001/rep5/m3/error22.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap001/rep5/input_ap001_rep5.D.tsv -o ./sim_output/ap001/rep5/m3
