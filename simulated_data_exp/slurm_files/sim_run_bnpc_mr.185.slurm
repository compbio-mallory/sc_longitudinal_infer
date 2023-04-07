#!/bin/bash
#SBATCH --job-name=bnpc185
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep3/m1/output185.out
#SBATCH --error=./sim_output/mc12/rep3/m1/error185.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep3/input_mc12_rep3.D.tsv -o ./sim_output/mc12/rep3/m1
