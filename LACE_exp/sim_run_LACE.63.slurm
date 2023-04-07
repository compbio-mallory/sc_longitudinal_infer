#!/bin/bash
#SBATCH --job-name=LACE63
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t4/rep4/output63.out
#SBATCH --error=sim_output/t4/rep4/error63.out
module load R
Rscript run_LACE.R ./sim_input/t4/rep4/t4_rep4.RData sim_output/t4/rep4/t4_rep4_inf.rdata sim_output/t4/rep4/t4_rep4_lgTree.xml
