#!/bin/bash
#SBATCH --job-name=LACE2
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_LACE_op/default/rep3/output2.out
#SBATCH --error=sim_LACE_op/default/rep3/error2.out
module load R
Rscript run_LACE.R ./sim_LACE_in/default/rep3/default_rep3.RData sim_LACE_op/default/rep3/inf_rep3.rdata sim_LACE_op/default/rep3/lgTree_rep3.xml
