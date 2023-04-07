#!/bin/bash
#SBATCH --job-name=LACE34
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/default/rep5/output34.out
#SBATCH --error=sim_output/default/rep5/error34.out
module load R
Rscript run_LACE.R ./sim_input/default/rep5/default_rep5.RData sim_output/default/rep5/default_rep5_inf.rdata sim_output/default/rep5/default_rep5_lgTree.xml
