#!/bin/bash
#SBATCH --job-name=LACE68
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t5/rep4/output68.out
#SBATCH --error=sim_output/t5/rep4/error68.out
module load R
Rscript run_LACE.R ./sim_input/t5/rep4/t5_rep4.RData sim_output/t5/rep4/t5_rep4_inf.rdata sim_output/t5/rep4/t5_rep4_lgTree.xml
