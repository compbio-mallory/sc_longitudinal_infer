#!/bin/bash
#SBATCH --job-name=LACE3
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_LACE_op/default/rep4/output3.out
#SBATCH --error=sim_LACE_op/default/rep4/error3.out
module load R
Rscript run_LACE.R ./sim_LACE_in/default/rep4/default_rep4.RData sim_LACE_op/default/rep4/inf_rep4.rdata sim_LACE_op/default/rep4/lgTree_rep4.xml
