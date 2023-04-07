#!/bin/bash
#SBATCH --job-name=LACE0
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_LACE_op/default/rep1/output0.out
#SBATCH --error=sim_LACE_op/default/rep1/error0.out
module load R
Rscript run_LACE.R ./sim_LACE_in/default/rep1/default_rep1.RData sim_LACE_op/default/rep1/inf_rep1.rdata sim_LACE_op/default/rep1/lgTree_rep1.xml
