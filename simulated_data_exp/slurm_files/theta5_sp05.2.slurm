#!/bin/bash
#SBATCH --job-name=bnpc2
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/theta5/rep3/m3/output2.out
#SBATCH --error=./sim_output/theta5/rep3/m3/error2.out
time python ../../BnpC/run_BnpC.py ./sim_input/theta5/rep3/input_theta5_rep3.D.tsv -o ./sim_output/theta5/rep3/m3
