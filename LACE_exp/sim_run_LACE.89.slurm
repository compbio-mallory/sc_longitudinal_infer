#!/bin/bash
#SBATCH --job-name=LACE89
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep5/output89.out
#SBATCH --error=sim_output/u1/rep5/error89.out
module load R
Rscript run_LACE.R ./sim_input/u1/rep5/u1_rep5.RData sim_output/u1/rep5/u1_rep5_inf.rdata sim_output/u1/rep5/u1_rep5_lgTree.xml
