#!/bin/bash
#SBATCH --job-name=bnpc130
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp4/rep2/m1/output130.out
#SBATCH --error=./sim_output/bp4/rep2/m1/error130.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp4/rep2/input_bp4_rep2.D.tsv -o ./sim_output/bp4/rep2/m1
