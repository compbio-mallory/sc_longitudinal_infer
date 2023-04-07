#!/bin/bash
#SBATCH --job-name=LACE42
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc2/rep3/output42.out
#SBATCH --error=sim_output/mc2/rep3/error42.out
module load R
Rscript run_LACE.R ./sim_input/mc2/rep3/mc2_rep3.RData sim_output/mc2/rep3/mc2_rep3_inf.rdata sim_output/mc2/rep3/mc2_rep3_lgTree.xml
