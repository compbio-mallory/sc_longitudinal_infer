#!/bin/bash
#SBATCH --job-name=LACE83
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u0/rep4/output83.out
#SBATCH --error=sim_output/u0/rep4/error83.out
module load R
Rscript run_LACE.R ./sim_input/u0/rep4/u0_rep4.RData sim_output/u0/rep4/u0_rep4_inf.rdata sim_output/u0/rep4/u0_rep4_lgTree.xml
