#!/bin/bash
#SBATCH --job-name=LACE87
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep3/output87.out
#SBATCH --error=sim_output/u1/rep3/error87.out
module load R
Rscript run_LACE.R ./sim_input/u1/rep3/u1_rep3.RData sim_output/u1/rep3/u1_rep3_inf.rdata sim_output/u1/rep3/u1_rep3_lgTree.xml
