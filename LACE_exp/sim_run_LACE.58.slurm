#!/bin/bash
#SBATCH --job-name=LACE58
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp5/rep4/output58.out
#SBATCH --error=sim_output/sp5/rep4/error58.out
module load R
Rscript run_LACE.R ./sim_input/sp5/rep4/sp5_rep4.RData sim_output/sp5/rep4/sp5_rep4_inf.rdata sim_output/sp5/rep4/sp5_rep4_lgTree.xml
