#!/bin/bash
#SBATCH --job-name=bnpc468
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep4/m4/output468.out
#SBATCH --error=./sim_output/u3/rep4/m4/error468.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep4/input_u3_rep4.D.tsv -o ./sim_output/u3/rep4/m4
