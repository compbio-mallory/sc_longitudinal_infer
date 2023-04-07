#!/bin/bash
#SBATCH --job-name=bnpc7
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/sp05/rep3/m3/output7.out
#SBATCH --error=./sim_output/sp05/rep3/m3/error7.out
time python ../../BnpC/run_BnpC.py ./sim_input/sp05/rep3/input_sp05_rep3.D.tsv -o ./sim_output/sp05/rep3/m3
