#!/bin/bash
#SBATCH --job-name=LACE49
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mp3/rep5/output49.out
#SBATCH --error=sim_output/mp3/rep5/error49.out
module load R
Rscript run_LACE.R ./sim_input/mp3/rep5/mp3_rep5.RData sim_output/mp3/rep5/mp3_rep5_inf.rdata sim_output/mp3/rep5/mp3_rep5_lgTree.xml
