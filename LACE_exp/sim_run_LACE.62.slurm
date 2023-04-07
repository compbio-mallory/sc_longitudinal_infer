#!/bin/bash
#SBATCH --job-name=LACE62
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/t4/rep3/output62.out
#SBATCH --error=sim_output/t4/rep3/error62.out
module load R
Rscript run_LACE.R ./sim_input/t4/rep3/t4_rep3.RData sim_output/t4/rep3/t4_rep3_inf.rdata sim_output/t4/rep3/t4_rep3_lgTree.xml
