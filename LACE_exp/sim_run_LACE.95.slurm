#!/bin/bash
#SBATCH --job-name=LACE95
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/v1v2/rep1/output95.out
#SBATCH --error=sim_output/v1v2/rep1/error95.out
module load R
Rscript run_LACE.R ./sim_input/v1v2/rep1/v1v2_rep1.RData sim_output/v1v2/rep1/v1v2_rep1_inf.rdata sim_output/v1v2/rep1/v1v2_rep1_lgTree.xml
