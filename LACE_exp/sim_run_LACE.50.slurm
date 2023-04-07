#!/bin/bash
#SBATCH --job-name=LACE50
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp05/rep1/output50.out
#SBATCH --error=sim_output/sp05/rep1/error50.out
module load R
Rscript run_LACE.R ./sim_input/sp05/rep1/sp05_rep1.RData sim_output/sp05/rep1/sp05_rep1_inf.rdata sim_output/sp05/rep1/sp05_rep1_lgTree.xml
