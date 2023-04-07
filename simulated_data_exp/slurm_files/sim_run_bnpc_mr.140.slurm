#!/bin/bash
#SBATCH --job-name=bnpc140
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp4/rep4/m1/output140.out
#SBATCH --error=./sim_output/bp4/rep4/m1/error140.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp4/rep4/input_bp4_rep4.D.tsv -o ./sim_output/bp4/rep4/m1
