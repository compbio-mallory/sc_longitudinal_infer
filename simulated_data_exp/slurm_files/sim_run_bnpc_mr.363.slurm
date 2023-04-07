#!/bin/bash
#SBATCH --job-name=bnpc363
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep3/m4/output363.out
#SBATCH --error=./sim_output/theta3/rep3/m4/error363.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep3/input_theta3_rep3.D.tsv -o ./sim_output/theta3/rep3/m4
