#!/bin/bash
#SBATCH --job-name=LACE60
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t4/rep1/output60.out
#SBATCH --error=sim_output/t4/rep1/error60.out
module load R
Rscript run_LACE.R ./sim_input/t4/rep1/t4_rep1.RData sim_output/t4/rep1/t4_rep1_inf.rdata sim_output/t4/rep1/t4_rep1_lgTree.xml
