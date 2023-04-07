#!/bin/bash
#SBATCH --job-name=LACE55
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp5/rep1/output55.out
#SBATCH --error=sim_output/sp5/rep1/error55.out
module load R
Rscript run_LACE.R ./sim_input/sp5/rep1/sp5_rep1.RData sim_output/sp5/rep1/sp5_rep1_inf.rdata sim_output/sp5/rep1/sp5_rep1_lgTree.xml
