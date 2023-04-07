#!/bin/bash
#SBATCH --job-name=LACE67
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t5/rep3/output67.out
#SBATCH --error=sim_output/t5/rep3/error67.out
module load R
Rscript run_LACE.R ./sim_input/t5/rep3/t5_rep3.RData sim_output/t5/rep3/t5_rep3_inf.rdata sim_output/t5/rep3/t5_rep3_lgTree.xml
