#!/bin/bash
#SBATCH --job-name=bnpc37
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap03/rep3/m3/output37.out
#SBATCH --error=./sim_output/ap03/rep3/m3/error37.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap03/rep3/input_ap03_rep3.D.tsv -o ./sim_output/ap03/rep3/m3
