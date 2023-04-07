#!/bin/bash
#SBATCH --job-name=LACE44
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/mc2/rep5/output44.out
#SBATCH --error=sim_output/mc2/rep5/error44.out
module load R
Rscript run_LACE.R ./sim_input/mc2/rep5/mc2_rep5.RData sim_output/mc2/rep5/mc2_rep5_inf.rdata sim_output/mc2/rep5/mc2_rep5_lgTree.xml
