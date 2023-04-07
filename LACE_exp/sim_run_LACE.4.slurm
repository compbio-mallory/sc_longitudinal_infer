#!/bin/bash
#SBATCH --job-name=LACE4
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap001/rep5/output4.out
#SBATCH --error=sim_output/ap001/rep5/error4.out
module load R
Rscript run_LACE.R ./sim_input/ap001/rep5/ap001_rep5.RData sim_output/ap001/rep5/ap001_rep5_inf.rdata sim_output/ap001/rep5/ap001_rep5_lgTree.xml
