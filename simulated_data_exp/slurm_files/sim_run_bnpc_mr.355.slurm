#!/bin/bash
#SBATCH --job-name=bnpc355
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta3/rep2/m1/output355.out
#SBATCH --error=./sim_output/theta3/rep2/m1/error355.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta3/rep2/input_theta3_rep2.D.tsv -o ./sim_output/theta3/rep2/m1
