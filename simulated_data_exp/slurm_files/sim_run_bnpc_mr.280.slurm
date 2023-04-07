#!/bin/bash
#SBATCH --job-name=bnpc280
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp5/rep2/m1/output280.out
#SBATCH --error=./sim_output/sp5/rep2/m1/error280.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp5/rep2/input_sp5_rep2.D.tsv -o ./sim_output/sp5/rep2/m1
