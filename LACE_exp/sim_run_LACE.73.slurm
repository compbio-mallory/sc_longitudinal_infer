#!/bin/bash
#SBATCH --job-name=LACE73
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta3/rep4/output73.out
#SBATCH --error=sim_output/theta3/rep4/error73.out
module load R
Rscript run_LACE.R ./sim_input/theta3/rep4/theta3_rep4.RData sim_output/theta3/rep4/theta3_rep4_inf.rdata sim_output/theta3/rep4/theta3_rep4_lgTree.xml
