#!/bin/bash
#SBATCH --job-name=LACE15
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp1/rep1/output15.out
#SBATCH --error=sim_output/bp1/rep1/error15.out
module load R
Rscript run_LACE.R ./sim_input/bp1/rep1/bp1_rep1.RData sim_output/bp1/rep1/bp1_rep1_inf.rdata sim_output/bp1/rep1/bp1_rep1_lgTree.xml
