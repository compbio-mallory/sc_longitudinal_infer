#!/bin/bash
#SBATCH --job-name=LACE39
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc12/rep5/output39.out
#SBATCH --error=sim_output/mc12/rep5/error39.out
module load R
Rscript run_LACE.R ./sim_input/mc12/rep5/mc12_rep5.RData sim_output/mc12/rep5/mc12_rep5_inf.rdata sim_output/mc12/rep5/mc12_rep5_lgTree.xml
