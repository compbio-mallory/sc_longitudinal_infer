#!/bin/bash
#SBATCH --job-name=bnpc474
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=./sim_output/u3/rep5/m5/output474.out
#SBATCH --error=./sim_output/u3/rep5/m5/error474.out
time python ../../BnpC/run_BnpC.py ./sim_input/u3/rep5/input_u3_rep5.D.tsv -o ./sim_output/u3/rep5/m5
