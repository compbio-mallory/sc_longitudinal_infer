#!/bin/bash
#SBATCH --job-name=LACE47
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mp3/rep3/output47.out
#SBATCH --error=sim_output/mp3/rep3/error47.out
module load R
Rscript run_LACE.R ./sim_input/mp3/rep3/mp3_rep3.RData sim_output/mp3/rep3/mp3_rep3_inf.rdata sim_output/mp3/rep3/mp3_rep3_lgTree.xml
