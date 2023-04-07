#!/bin/bash
#SBATCH --job-name=LACE31
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/default/rep2/output31.out
#SBATCH --error=sim_output/default/rep2/error31.out
module load R
Rscript run_LACE.R ./sim_input/default/rep2/default_rep2.RData sim_output/default/rep2/default_rep2_inf.rdata sim_output/default/rep2/default_rep2_lgTree.xml
