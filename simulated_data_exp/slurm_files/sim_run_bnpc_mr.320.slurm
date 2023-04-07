#!/bin/bash
#SBATCH --job-name=bnpc320
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/t4/rep5/m1/output320.out
#SBATCH --error=./sim_output/t4/rep5/m1/error320.out
time python ../../BnpC/run_BnpC.py ./sim_input/t4/rep5/t4_rep5.D.tsv -o ./sim_output/t4/rep5/m1
