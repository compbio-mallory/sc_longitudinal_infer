#!/bin/bash
#SBATCH --job-name=LACE17
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp1/rep3/output17.out
#SBATCH --error=sim_output/bp1/rep3/error17.out
module load R
Rscript run_LACE.R ./sim_input/bp1/rep3/bp1_rep3.RData sim_output/bp1/rep3/bp1_rep3_inf.rdata sim_output/bp1/rep3/bp1_rep3_lgTree.xml
