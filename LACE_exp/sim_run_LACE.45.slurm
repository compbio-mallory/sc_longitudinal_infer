#!/bin/bash
#SBATCH --job-name=LACE45
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mp3/rep1/output45.out
#SBATCH --error=sim_output/mp3/rep1/error45.out
module load R
Rscript run_LACE.R ./sim_input/mp3/rep1/mp3_rep1.RData sim_output/mp3/rep1/mp3_rep1_inf.rdata sim_output/mp3/rep1/mp3_rep1_lgTree.xml
