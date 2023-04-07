#!/bin/bash
#SBATCH --job-name=bnpc357
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep2/m3/output357.out
#SBATCH --error=./sim_output/theta3/rep2/m3/error357.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep2/input_theta3_rep2.D.tsv -o ./sim_output/theta3/rep2/m3
