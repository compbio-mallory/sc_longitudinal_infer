#!/bin/bash
#SBATCH --job-name=bnpc335
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t5/rep3/m1/output335.out
#SBATCH --error=./sim_output/t5/rep3/m1/error335.out
time python ../../BnpC/run_BnpC.py ./sim_input/t5/rep3/t5_rep3.D.tsv -o ./sim_output/t5/rep3/m1
