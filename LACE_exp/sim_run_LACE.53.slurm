#!/bin/bash
#SBATCH --job-name=LACE53
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp05/rep4/output53.out
#SBATCH --error=sim_output/sp05/rep4/error53.out
module load R
Rscript run_LACE.R ./sim_input/sp05/rep4/sp05_rep4.RData sim_output/sp05/rep4/sp05_rep4_inf.rdata sim_output/sp05/rep4/sp05_rep4_lgTree.xml
