#!/bin/bash
#SBATCH --job-name=LACE36
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc12/rep2/output36.out
#SBATCH --error=sim_output/mc12/rep2/error36.out
module load R
Rscript run_LACE.R ./sim_input/mc12/rep2/mc12_rep2.RData sim_output/mc12/rep2/mc12_rep2_inf.rdata sim_output/mc12/rep2/mc12_rep2_lgTree.xml
