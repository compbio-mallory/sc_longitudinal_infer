#!/bin/bash
#SBATCH --job-name=LACE99
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/v1v2/rep5/output99.out
#SBATCH --error=sim_output/v1v2/rep5/error99.out
module load R
Rscript run_LACE.R ./sim_input/v1v2/rep5/v1v2_rep5.RData sim_output/v1v2/rep5/v1v2_rep5_inf.rdata sim_output/v1v2/rep5/v1v2_rep5_lgTree.xml
