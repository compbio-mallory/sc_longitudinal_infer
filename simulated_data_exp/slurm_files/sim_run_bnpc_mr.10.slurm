#!/bin/bash
#SBATCH --job-name=bnpc10
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap001/rep3/m1/output10.out
#SBATCH --error=./sim_output/ap001/rep3/m1/error10.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap001/rep3/input_ap001_rep3.D.tsv -o ./sim_output/ap001/rep3/m1
