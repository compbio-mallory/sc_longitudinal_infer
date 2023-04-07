#!/bin/bash
#SBATCH --job-name=LACE23
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp3/rep4/output23.out
#SBATCH --error=sim_output/bp3/rep4/error23.out
module load R
Rscript run_LACE.R ./sim_input/bp3/rep4/bp3_rep4.RData sim_output/bp3/rep4/bp3_rep4_inf.rdata sim_output/bp3/rep4/bp3_rep4_lgTree.xml
