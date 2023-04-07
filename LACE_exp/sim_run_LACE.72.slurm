#!/bin/bash
#SBATCH --job-name=LACE72
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta3/rep3/output72.out
#SBATCH --error=sim_output/theta3/rep3/error72.out
module load R
Rscript run_LACE.R ./sim_input/theta3/rep3/theta3_rep3.RData sim_output/theta3/rep3/theta3_rep3_inf.rdata sim_output/theta3/rep3/theta3_rep3_lgTree.xml
