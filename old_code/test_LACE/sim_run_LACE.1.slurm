#!/bin/bash
#SBATCH --job-name=LACE1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_LACE_op/default/rep2/output1.out
#SBATCH --error=sim_LACE_op/default/rep2/error1.out
module load R
Rscript run_LACE.R ./sim_LACE_in/default/rep2/default_rep2.RData sim_LACE_op/default/rep2/inf_rep2.rdata sim_LACE_op/default/rep2/lgTree_rep2.xml
