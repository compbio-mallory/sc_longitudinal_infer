#!/bin/bash
#SBATCH --job-name=bnpc125
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp4/rep1/m1/output125.out
#SBATCH --error=./sim_output/bp4/rep1/m1/error125.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp4/rep1/input_bp4_rep1.D.tsv -o ./sim_output/bp4/rep1/m1
