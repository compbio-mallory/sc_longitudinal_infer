#!/bin/bash
#SBATCH --job-name=bnpc181
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep2/m2/output181.out
#SBATCH --error=./sim_output/mc12/rep2/m2/error181.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep2/input_mc12_rep2.D.tsv -o ./sim_output/mc12/rep2/m2
