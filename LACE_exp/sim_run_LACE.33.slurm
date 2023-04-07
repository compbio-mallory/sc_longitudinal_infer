#!/bin/bash
#SBATCH --job-name=LACE33
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/default/rep4/output33.out
#SBATCH --error=sim_output/default/rep4/error33.out
module load R
Rscript run_LACE.R ./sim_input/default/rep4/default_rep4.RData sim_output/default/rep4/default_rep4_inf.rdata sim_output/default/rep4/default_rep4_lgTree.xml
