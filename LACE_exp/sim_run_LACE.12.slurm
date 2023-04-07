#!/bin/bash
#SBATCH --job-name=LACE12
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap05/rep3/output12.out
#SBATCH --error=sim_output/ap05/rep3/error12.out
module load R
Rscript run_LACE.R ./sim_input/ap05/rep3/ap05_rep3.RData sim_output/ap05/rep3/ap05_rep3_inf.rdata sim_output/ap05/rep3/ap05_rep3_lgTree.xml
