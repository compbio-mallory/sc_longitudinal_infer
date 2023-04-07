#!/bin/bash
#SBATCH --job-name=LACE10
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap05/rep1/output10.out
#SBATCH --error=sim_output/ap05/rep1/error10.out
module load R
Rscript run_LACE.R ./sim_input/ap05/rep1/ap05_rep1.RData sim_output/ap05/rep1/ap05_rep1_inf.rdata sim_output/ap05/rep1/ap05_rep1_lgTree.xml
