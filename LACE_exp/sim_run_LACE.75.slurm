#!/bin/bash
#SBATCH --job-name=LACE75
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta5/rep1/output75.out
#SBATCH --error=sim_output/theta5/rep1/error75.out
module load R
Rscript run_LACE.R ./sim_input/theta5/rep1/theta5_rep1.RData sim_output/theta5/rep1/theta5_rep1_inf.rdata sim_output/theta5/rep1/theta5_rep1_lgTree.xml
