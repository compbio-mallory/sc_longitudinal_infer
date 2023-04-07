#!/bin/bash
#SBATCH --job-name=LACE51
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp05/rep2/output51.out
#SBATCH --error=sim_output/sp05/rep2/error51.out
module load R
Rscript run_LACE.R ./sim_input/sp05/rep2/sp05_rep2.RData sim_output/sp05/rep2/sp05_rep2_inf.rdata sim_output/sp05/rep2/sp05_rep2_lgTree.xml
