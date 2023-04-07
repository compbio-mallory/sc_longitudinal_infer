#!/bin/bash
#SBATCH --job-name=LACE9
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap03/rep5/output9.out
#SBATCH --error=sim_output/ap03/rep5/error9.out
module load R
Rscript run_LACE.R ./sim_input/ap03/rep5/ap03_rep5.RData sim_output/ap03/rep5/ap03_rep5_inf.rdata sim_output/ap03/rep5/ap03_rep5_lgTree.xml
