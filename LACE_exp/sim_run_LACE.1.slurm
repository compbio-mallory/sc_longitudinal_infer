#!/bin/bash
#SBATCH --job-name=LACE1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/ap001/rep2/output1.out
#SBATCH --error=sim_output/ap001/rep2/error1.out
module load R
Rscript run_LACE.R ./sim_input/ap001/rep2/ap001_rep2.RData sim_output/ap001/rep2/ap001_rep2_inf.rdata sim_output/ap001/rep2/ap001_rep2_lgTree.xml
