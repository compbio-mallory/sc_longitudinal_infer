#!/bin/bash
#SBATCH --job-name=LACE57
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp5/rep3/output57.out
#SBATCH --error=sim_output/sp5/rep3/error57.out
module load R
Rscript run_LACE.R ./sim_input/sp5/rep3/sp5_rep3.RData sim_output/sp5/rep3/sp5_rep3_inf.rdata sim_output/sp5/rep3/sp5_rep3_lgTree.xml
