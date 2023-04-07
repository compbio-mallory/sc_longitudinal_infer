#!/bin/bash
#SBATCH --job-name=LACE32
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/default/rep3/output32.out
#SBATCH --error=sim_output/default/rep3/error32.out
module load R
Rscript run_LACE.R ./sim_input/default/rep3/default_rep3.RData sim_output/default/rep3/default_rep3_inf.rdata sim_output/default/rep3/default_rep3_lgTree.xml
