#!/bin/bash
#SBATCH --job-name=bnpc180
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc12/rep2/m1/output180.out
#SBATCH --error=./sim_output/mc12/rep2/m1/error180.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc12/rep2/input_mc12_rep2.D.tsv -o ./sim_output/mc12/rep2/m1
