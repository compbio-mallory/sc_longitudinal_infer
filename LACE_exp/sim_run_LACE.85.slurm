#!/bin/bash
#SBATCH --job-name=LACE85
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep1/output85.out
#SBATCH --error=sim_output/u1/rep1/error85.out
module load R
Rscript run_LACE.R ./sim_input/u1/rep1/u1_rep1.RData sim_output/u1/rep1/u1_rep1_inf.rdata sim_output/u1/rep1/u1_rep1_lgTree.xml
