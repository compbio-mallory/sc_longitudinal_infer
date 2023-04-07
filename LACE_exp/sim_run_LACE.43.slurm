#!/bin/bash
#SBATCH --job-name=LACE43
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc2/rep4/output43.out
#SBATCH --error=sim_output/mc2/rep4/error43.out
module load R
Rscript run_LACE.R ./sim_input/mc2/rep4/mc2_rep4.RData sim_output/mc2/rep4/mc2_rep4_inf.rdata sim_output/mc2/rep4/mc2_rep4_lgTree.xml
