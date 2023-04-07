#!/bin/bash
#SBATCH --job-name=LACE78
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta5/rep4/output78.out
#SBATCH --error=sim_output/theta5/rep4/error78.out
module load R
Rscript run_LACE.R ./sim_input/theta5/rep4/theta5_rep4.RData sim_output/theta5/rep4/theta5_rep4_inf.rdata sim_output/theta5/rep4/theta5_rep4_lgTree.xml
