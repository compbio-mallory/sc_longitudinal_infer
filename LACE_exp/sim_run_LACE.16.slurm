#!/bin/bash
#SBATCH --job-name=LACE16
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/bp1/rep2/output16.out
#SBATCH --error=sim_output/bp1/rep2/error16.out
module load R
Rscript run_LACE.R ./sim_input/bp1/rep2/bp1_rep2.RData sim_output/bp1/rep2/bp1_rep2_inf.rdata sim_output/bp1/rep2/bp1_rep2_lgTree.xml
