#!/bin/bash
#SBATCH --job-name=bnpc398
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep5/m4/output398.out
#SBATCH --error=./sim_output/theta5/rep5/m4/error398.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep5/input_theta5_rep5.D.tsv -o ./sim_output/theta5/rep5/m4
