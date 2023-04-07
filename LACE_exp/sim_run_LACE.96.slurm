#!/bin/bash
#SBATCH --job-name=LACE96
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/v1v2/rep2/output96.out
#SBATCH --error=sim_output/v1v2/rep2/error96.out
module load R
Rscript run_LACE.R ./sim_input/v1v2/rep2/v1v2_rep2.RData sim_output/v1v2/rep2/v1v2_rep2_inf.rdata sim_output/v1v2/rep2/v1v2_rep2_lgTree.xml
