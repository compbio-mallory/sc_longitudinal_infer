#!/bin/bash
#SBATCH --job-name=LACE98
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/v1v2/rep4/output98.out
#SBATCH --error=sim_output/v1v2/rep4/error98.out
module load R
Rscript run_LACE.R ./sim_input/v1v2/rep4/v1v2_rep4.RData sim_output/v1v2/rep4/v1v2_rep4_inf.rdata sim_output/v1v2/rep4/v1v2_rep4_lgTree.xml
