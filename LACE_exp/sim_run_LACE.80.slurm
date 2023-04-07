#!/bin/bash
#SBATCH --job-name=LACE80
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u0/rep1/output80.out
#SBATCH --error=sim_output/u0/rep1/error80.out
module load R
Rscript run_LACE.R ./sim_input/u0/rep1/u0_rep1.RData sim_output/u0/rep1/u0_rep1_inf.rdata sim_output/u0/rep1/u0_rep1_lgTree.xml
