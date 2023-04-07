#!/bin/bash
#SBATCH --job-name=LACE4
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_LACE_op/default/rep5/output4.out
#SBATCH --error=sim_LACE_op/default/rep5/error4.out
module load R
Rscript run_LACE.R ./sim_LACE_in/default/rep5/default_rep5.RData sim_LACE_op/default/rep5/inf_rep5.rdata sim_LACE_op/default/rep5/lgTree_rep5.xml
