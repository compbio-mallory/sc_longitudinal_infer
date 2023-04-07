#!/bin/bash
#SBATCH --job-name=LACE11
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap05/rep2/output11.out
#SBATCH --error=sim_output/ap05/rep2/error11.out
module load R
Rscript run_LACE.R ./sim_input/ap05/rep2/ap05_rep2.RData sim_output/ap05/rep2/ap05_rep2_inf.rdata sim_output/ap05/rep2/ap05_rep2_lgTree.xml
