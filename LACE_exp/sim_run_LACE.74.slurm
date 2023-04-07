#!/bin/bash
#SBATCH --job-name=LACE74
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta3/rep5/output74.out
#SBATCH --error=sim_output/theta3/rep5/error74.out
module load R
Rscript run_LACE.R ./sim_input/theta3/rep5/theta3_rep5.RData sim_output/theta3/rep5/theta3_rep5_inf.rdata sim_output/theta3/rep5/theta3_rep5_lgTree.xml
