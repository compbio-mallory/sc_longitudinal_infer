#!/bin/bash
#SBATCH --job-name=LACE76
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta5/rep2/output76.out
#SBATCH --error=sim_output/theta5/rep2/error76.out
module load R
Rscript run_LACE.R ./sim_input/theta5/rep2/theta5_rep2.RData sim_output/theta5/rep2/theta5_rep2_inf.rdata sim_output/theta5/rep2/theta5_rep2_lgTree.xml
