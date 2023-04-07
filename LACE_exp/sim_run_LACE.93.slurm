#!/bin/bash
#SBATCH --job-name=LACE93
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/u3/rep4/output93.out
#SBATCH --error=sim_output/u3/rep4/error93.out
module load R
Rscript run_LACE.R ./sim_input/u3/rep4/u3_rep4.RData sim_output/u3/rep4/u3_rep4_inf.rdata sim_output/u3/rep4/u3_rep4_lgTree.xml
