#!/bin/bash
#SBATCH --job-name=LACE46
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mp3/rep2/output46.out
#SBATCH --error=sim_output/mp3/rep2/error46.out
module load R
Rscript run_LACE.R ./sim_input/mp3/rep2/mp3_rep2.RData sim_output/mp3/rep2/mp3_rep2_inf.rdata sim_output/mp3/rep2/mp3_rep2_lgTree.xml
