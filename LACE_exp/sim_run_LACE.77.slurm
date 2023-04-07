#!/bin/bash
#SBATCH --job-name=LACE77
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/theta5/rep3/output77.out
#SBATCH --error=sim_output/theta5/rep3/error77.out
module load R
Rscript run_LACE.R ./sim_input/theta5/rep3/theta5_rep3.RData sim_output/theta5/rep3/theta5_rep3_inf.rdata sim_output/theta5/rep3/theta5_rep3_lgTree.xml
