#!/bin/bash
#SBATCH --job-name=LACE41
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc2/rep2/output41.out
#SBATCH --error=sim_output/mc2/rep2/error41.out
module load R
Rscript run_LACE.R ./sim_input/mc2/rep2/mc2_rep2.RData sim_output/mc2/rep2/mc2_rep2_inf.rdata sim_output/mc2/rep2/mc2_rep2_lgTree.xml
