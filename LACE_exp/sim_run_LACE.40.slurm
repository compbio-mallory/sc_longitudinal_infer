#!/bin/bash
#SBATCH --job-name=LACE40
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc2/rep1/output40.out
#SBATCH --error=sim_output/mc2/rep1/error40.out
module load R
Rscript run_LACE.R ./sim_input/mc2/rep1/mc2_rep1.RData sim_output/mc2/rep1/mc2_rep1_inf.rdata sim_output/mc2/rep1/mc2_rep1_lgTree.xml
