#!/bin/bash
#SBATCH --job-name=LACE82
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u0/rep3/output82.out
#SBATCH --error=sim_output/u0/rep3/error82.out
module load R
Rscript run_LACE.R ./sim_input/u0/rep3/u0_rep3.RData sim_output/u0/rep3/u0_rep3_inf.rdata sim_output/u0/rep3/u0_rep3_lgTree.xml
