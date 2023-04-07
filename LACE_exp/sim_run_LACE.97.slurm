#!/bin/bash
#SBATCH --job-name=LACE97
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/v1v2/rep3/output97.out
#SBATCH --error=sim_output/v1v2/rep3/error97.out
module load R
Rscript run_LACE.R ./sim_input/v1v2/rep3/v1v2_rep3.RData sim_output/v1v2/rep3/v1v2_rep3_inf.rdata sim_output/v1v2/rep3/v1v2_rep3_lgTree.xml
