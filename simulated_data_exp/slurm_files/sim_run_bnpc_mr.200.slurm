#!/bin/bash
#SBATCH --job-name=bnpc200
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/mc2/rep1/m1/output200.out
#SBATCH --error=./sim_output/mc2/rep1/m1/error200.out
time python ../../BnpC/run_BnpC.py ./sim_input/mc2/rep1/input_mc2_rep1.D.tsv -o ./sim_output/mc2/rep1/m1
