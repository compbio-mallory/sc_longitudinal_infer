#!/bin/bash
#SBATCH --job-name=LACE48
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mp3/rep4/output48.out
#SBATCH --error=sim_output/mp3/rep4/error48.out
module load R
Rscript run_LACE.R ./sim_input/mp3/rep4/mp3_rep4.RData sim_output/mp3/rep4/mp3_rep4_inf.rdata sim_output/mp3/rep4/mp3_rep4_lgTree.xml
