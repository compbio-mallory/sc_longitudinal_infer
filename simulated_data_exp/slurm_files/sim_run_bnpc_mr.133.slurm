#!/bin/bash
#SBATCH --job-name=bnpc133
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/bp4/rep2/m4/output133.out
#SBATCH --error=./sim_output/bp4/rep2/m4/error133.out
time python ../../BnpC/run_BnpC.py ./sim_input/bp4/rep2/input_bp4_rep2.D.tsv -o ./sim_output/bp4/rep2/m4
