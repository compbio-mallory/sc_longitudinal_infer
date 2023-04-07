#!/bin/bash
#SBATCH --job-name=bnpc0
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=4GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=bnpc_output/m5/output0.out
#SBATCH --error=bnpc_output/m5/error0.out
time python ../../../BnpC/run_BnpC.py SA501/bnpcInput.D.tsv -o bnpc_output/m5/
