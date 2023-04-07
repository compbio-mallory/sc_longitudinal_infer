#!/bin/bash
#SBATCH --job-name=LACE65
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t5/rep1/output65.out
#SBATCH --error=sim_output/t5/rep1/error65.out
module load R
Rscript run_LACE.R ./sim_input/t5/rep1/t5_rep1.RData sim_output/t5/rep1/t5_rep1_inf.rdata sim_output/t5/rep1/t5_rep1_lgTree.xml
