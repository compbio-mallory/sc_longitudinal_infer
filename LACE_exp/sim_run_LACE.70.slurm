#!/bin/bash
#SBATCH --job-name=LACE70
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta3/rep1/output70.out
#SBATCH --error=sim_output/theta3/rep1/error70.out
module load R
Rscript run_LACE.R ./sim_input/theta3/rep1/theta3_rep1.RData sim_output/theta3/rep1/theta3_rep1_inf.rdata sim_output/theta3/rep1/theta3_rep1_lgTree.xml
