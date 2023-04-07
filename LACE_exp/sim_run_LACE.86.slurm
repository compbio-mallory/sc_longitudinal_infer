#!/bin/bash
#SBATCH --job-name=LACE86
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u1/rep2/output86.out
#SBATCH --error=sim_output/u1/rep2/error86.out
module load R
Rscript run_LACE.R ./sim_input/u1/rep2/u1_rep2.RData sim_output/u1/rep2/u1_rep2_inf.rdata sim_output/u1/rep2/u1_rep2_lgTree.xml
