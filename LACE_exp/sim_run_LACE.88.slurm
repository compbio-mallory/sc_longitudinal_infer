#!/bin/bash
#SBATCH --job-name=LACE88
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep4/output88.out
#SBATCH --error=sim_output/u1/rep4/error88.out
module load R
Rscript run_LACE.R ./sim_input/u1/rep4/u1_rep4.RData sim_output/u1/rep4/u1_rep4_inf.rdata sim_output/u1/rep4/u1_rep4_lgTree.xml
