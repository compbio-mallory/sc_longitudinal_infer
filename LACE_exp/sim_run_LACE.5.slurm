#!/bin/bash
#SBATCH --job-name=LACE5
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap03/rep1/output5.out
#SBATCH --error=sim_output/ap03/rep1/error5.out
module load R
Rscript run_LACE.R ./sim_input/ap03/rep1/ap03_rep1.RData sim_output/ap03/rep1/ap03_rep1_inf.rdata sim_output/ap03/rep1/ap03_rep1_lgTree.xml
