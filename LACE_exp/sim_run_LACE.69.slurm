#!/bin/bash
#SBATCH --job-name=LACE69
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t5/rep5/output69.out
#SBATCH --error=sim_output/t5/rep5/error69.out
module load R
Rscript run_LACE.R ./sim_input/t5/rep5/t5_rep5.RData sim_output/t5/rep5/t5_rep5_inf.rdata sim_output/t5/rep5/t5_rep5_lgTree.xml
