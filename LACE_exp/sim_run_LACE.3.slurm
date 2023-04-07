#!/bin/bash
#SBATCH --job-name=LACE3
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap001/rep4/output3.out
#SBATCH --error=sim_output/ap001/rep4/error3.out
module load R
Rscript run_LACE.R ./sim_input/ap001/rep4/ap001_rep4.RData sim_output/ap001/rep4/ap001_rep4_inf.rdata sim_output/ap001/rep4/ap001_rep4_lgTree.xml
