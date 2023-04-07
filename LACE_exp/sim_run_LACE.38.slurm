#!/bin/bash
#SBATCH --job-name=LACE38
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc12/rep4/output38.out
#SBATCH --error=sim_output/mc12/rep4/error38.out
module load R
Rscript run_LACE.R ./sim_input/mc12/rep4/mc12_rep4.RData sim_output/mc12/rep4/mc12_rep4_inf.rdata sim_output/mc12/rep4/mc12_rep4_lgTree.xml
