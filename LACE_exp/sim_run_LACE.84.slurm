#!/bin/bash
#SBATCH --job-name=LACE84
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u0/rep5/output84.out
#SBATCH --error=sim_output/u0/rep5/error84.out
module load R
Rscript run_LACE.R ./sim_input/u0/rep5/u0_rep5.RData sim_output/u0/rep5/u0_rep5_inf.rdata sim_output/u0/rep5/u0_rep5_lgTree.xml
