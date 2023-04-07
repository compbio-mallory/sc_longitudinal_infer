#!/bin/bash
#SBATCH --job-name=LACE56
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=8GB
#SBATCH -p fan_q
#SBATCH --mail-type=ALL
#SBATCH --output=sim_output/sp5/rep2/output56.out
#SBATCH --error=sim_output/sp5/rep2/error56.out
module load R
Rscript run_LACE.R ./sim_input/sp5/rep2/sp5_rep2.RData sim_output/sp5/rep2/sp5_rep2_inf.rdata sim_output/sp5/rep2/sp5_rep2_lgTree.xml
