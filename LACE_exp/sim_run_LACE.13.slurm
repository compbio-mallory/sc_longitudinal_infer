#!/bin/bash
#SBATCH --job-name=LACE13
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap05/rep4/output13.out
#SBATCH --error=sim_output/ap05/rep4/error13.out
module load R
Rscript run_LACE.R ./sim_input/ap05/rep4/ap05_rep4.RData sim_output/ap05/rep4/ap05_rep4_inf.rdata sim_output/ap05/rep4/ap05_rep4_lgTree.xml
