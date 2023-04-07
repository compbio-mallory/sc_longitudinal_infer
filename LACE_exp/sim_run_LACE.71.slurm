#!/bin/bash
#SBATCH --job-name=LACE71
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta3/rep2/output71.out
#SBATCH --error=sim_output/theta3/rep2/error71.out
module load R
Rscript run_LACE.R ./sim_input/theta3/rep2/theta3_rep2.RData sim_output/theta3/rep2/theta3_rep2_inf.rdata sim_output/theta3/rep2/theta3_rep2_lgTree.xml
