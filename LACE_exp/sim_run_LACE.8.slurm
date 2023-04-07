#!/bin/bash
#SBATCH --job-name=LACE8
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap03/rep4/output8.out
#SBATCH --error=sim_output/ap03/rep4/error8.out
module load R
Rscript run_LACE.R ./sim_input/ap03/rep4/ap03_rep4.RData sim_output/ap03/rep4/ap03_rep4_inf.rdata sim_output/ap03/rep4/ap03_rep4_lgTree.xml
