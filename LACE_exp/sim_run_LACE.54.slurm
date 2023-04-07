#!/bin/bash
#SBATCH --job-name=LACE54
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp05/rep5/output54.out
#SBATCH --error=sim_output/sp05/rep5/error54.out
module load R
Rscript run_LACE.R ./sim_input/sp05/rep5/sp05_rep5.RData sim_output/sp05/rep5/sp05_rep5_inf.rdata sim_output/sp05/rep5/sp05_rep5_lgTree.xml
