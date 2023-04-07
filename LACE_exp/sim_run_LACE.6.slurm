#!/bin/bash
#SBATCH --job-name=LACE6
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap03/rep2/output6.out
#SBATCH --error=sim_output/ap03/rep2/error6.out
module load R
Rscript run_LACE.R ./sim_input/ap03/rep2/ap03_rep2.RData sim_output/ap03/rep2/ap03_rep2_inf.rdata sim_output/ap03/rep2/ap03_rep2_lgTree.xml
