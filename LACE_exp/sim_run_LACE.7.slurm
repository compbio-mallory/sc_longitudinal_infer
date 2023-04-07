#!/bin/bash
#SBATCH --job-name=LACE7
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap03/rep3/output7.out
#SBATCH --error=sim_output/ap03/rep3/error7.out
module load R
Rscript run_LACE.R ./sim_input/ap03/rep3/ap03_rep3.RData sim_output/ap03/rep3/ap03_rep3_inf.rdata sim_output/ap03/rep3/ap03_rep3_lgTree.xml
