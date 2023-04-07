#!/bin/bash
#SBATCH --job-name=bnpc1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/ap001/rep1/m2/output1.out
#SBATCH --error=./sim_output/ap001/rep1/m2/error1.out
time python ../../BnpC/run_BnpC.py ./sim_input/ap001/rep1/input_ap001_rep1.D.tsv -o ./sim_output/ap001/rep1/m2
