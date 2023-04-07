#!/bin/bash
#SBATCH --job-name=LACE30
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/default/rep1/output30.out
#SBATCH --error=sim_output/default/rep1/error30.out
module load R
Rscript run_LACE.R ./sim_input/default/rep1/default_rep1.RData sim_output/default/rep1/default_rep1_inf.rdata sim_output/default/rep1/default_rep1_lgTree.xml
