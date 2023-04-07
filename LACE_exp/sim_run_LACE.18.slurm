#!/bin/bash
#SBATCH --job-name=LACE18
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp1/rep4/output18.out
#SBATCH --error=sim_output/bp1/rep4/error18.out
module load R
Rscript run_LACE.R ./sim_input/bp1/rep4/bp1_rep4.RData sim_output/bp1/rep4/bp1_rep4_inf.rdata sim_output/bp1/rep4/bp1_rep4_lgTree.xml
