#!/bin/bash
#SBATCH --job-name=LACE35
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc12/rep1/output35.out
#SBATCH --error=sim_output/mc12/rep1/error35.out
module load R
Rscript run_LACE.R ./sim_input/mc12/rep1/mc12_rep1.RData sim_output/mc12/rep1/mc12_rep1_inf.rdata sim_output/mc12/rep1/mc12_rep1_lgTree.xml
